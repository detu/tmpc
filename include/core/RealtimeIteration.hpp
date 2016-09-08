#pragma once

#include "Trajectory.hpp"
#include "matrix.hpp"
#include "../qp/qp.hpp"

#include <stdexcept>
//#include <type_traits>

namespace tmpc
{
	template <typename OCP,
		template <unsigned, unsigned, unsigned, unsigned> class QPSolver_>
	class RealtimeIteration
	{
		static auto const NX = OCP::NX;
		static auto const NU = OCP::NU;
		static auto const NC = OCP::NC;
		static auto const NCT = OCP::NCT;

	public:
		typedef QPSolver_<NX, NU, NC, NCT> QPSolver;
		typedef Trajectory<NX, NU> WorkingPoint;

		RealtimeIteration(OCP const& ocp, QPSolver& solver, WorkingPoint const& working_point)
		:	_ocp(ocp)
		,	_QP(working_point.nT())
		,	workingPoint_(working_point)
		,	solution_(working_point.nT())
		,	solver_(solver)
		,	_prepared(true)
		{
			InitQP();
		}

		RealtimeIteration(RealtimeIteration const&) = delete;

		/**
		 * \brief Number of intervals.
		 */
		std::size_t nT() const { return workingPoint_.nT(); }

		// Feed current state x0, get back control input u.
		template <typename StateVector>
		decltype(auto) Feedback(const StateVector& x0)
		{
			if (!_prepared)
				throw std::logic_error("ModelPredictiveController::Feedback(): controller is not prepared.");

			/** embed current initial value */
			auto const w0 = eval(x0 - workingPoint_.get_x(0));
			_QP.set_x_min(0, w0);
			_QP.set_x_max(0, w0);

			/** solve QP */
			solver_.Solve(_QP, solution_);

			// Add QP step to the working point.
			{
				for (std::size_t i = 0; i < nT() + 1; ++i)
					workingPoint_.set_x(i, workingPoint_.get_x(i) + solution_.get_x(i));

				for (std::size_t i = 0; i < nT(); ++i)
					workingPoint_.set_u(i, workingPoint_.get_u(i) + solution_.get_u(i));
			}

			_prepared = false;

			// Return the calculated control input.
			return workingPoint_.get_u(0);
		}

		void Preparation()
		{
			if (_prepared)
				throw std::logic_error("ModelPredictiveController::Preparation(): controller is already prepared.");

			shift(workingPoint_);

			// Calculate new QP.
			UpdateQP();

			_prepared = true;
		}

		WorkingPoint const& getWorkingPoint() const { return workingPoint_; }

		double getLevenbergMarquardt() const { return _levenbergMarquardt; }
		void setLevenbergMarquardt(double val) { _levenbergMarquardt = val; }

		/*
		unsigned nT() const { return _ocp.getNumberOfIntervals(); }
		unsigned nU() const	noexcept { return _Nu; }
		unsigned nX() const	noexcept { return _Nx; }
		unsigned nZ() const noexcept { return _Nz; }
		*/

	private:
		typedef typename QPSolver::Solution Solution;
		typedef typename QPSolver::Problem QP;

		/**
		 * \brief Interface for setting MPC stage data.
		 * OCP::InitStage() and OCP::UpdateStage() use this interface to set Hessian of the stage cost,
		 * gradient of the stage cost, linearized stage constraints, state and input bounds.
		 *
		 * Hessian and gradient of Lagrange term:
		 * H = [Q, S
		 *      S, R]
		 * g = [q
		 *      r]
		 *
		 * NOTE:
		 * Potential issue with this implementation: it is easy to not to set one of the properties.
		 * But a single function with 11 arguments would look ugly.
		 * On the other hand, re-evaluation of the stage may be not necessary for some problems,
		 * and previous values would be kept. This allows more efficient implementations in some cases.
		 */
		class Stage
		{
		public:
			Stage(std::size_t i, WorkingPoint& wp, QP& qp, double lm)
			:	i_(i), workingPoint_(wp), qp_(qp), levenbergMarquardt_(lm) {
			}

			Stage(Stage const&) = delete;

			/**
			 * \brief Returns stage index.
			 */
			std::size_t get_index() const { return i_; }

			/**
			 * \brief Returns stage state.
			 */
			decltype(auto) get_x() const { return workingPoint_.get_x(i_); }

			/**
			 * \brief Returns stage input.
			 */
			decltype(auto) get_u() const { return workingPoint_.get_u(i_); }

			template <typename Matrix>
			void set_Q(Matrix const &Q)	{
				qp_.set_Q(i_, Q + levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_R(Matrix const &R) {
				qp_.set_R(i_, R + levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_S(Matrix const &S) { qp_.set_S(i_, S); }

			template <typename Vector>
			void set_q(Vector const &q) { qp_.set_q(i_, q); }

			template <typename Vector>
			void set_r(Vector const &r) { qp_.set_r(i_, r); }

			template <typename Matrix>
			void set_A(Matrix const &A) { qp_.set_A(i_, A); }

			template <typename Matrix>
			void set_B(Matrix const &B) { qp_.set_B(i_, B); }

			template <typename Vector>
			void set_x_next(Vector const& x_next)
			{
				// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
				// c = f(z_k) - x_{k+1}
				qp_.set_b(i_, x_next - workingPoint_.get_x(i_ + 1));
			}

			template <typename Matrix>
			void set_C(Matrix const &C) { qp_.set_C(i_, C); }

			template <typename Matrix>
			void set_D(Matrix const &D) { qp_.set_D(i_, D); }

			template <typename Vector>
			void set_d_min(Vector const& d_min) { qp_.set_d_min(i_, d_min); }

			template <typename Vector>
			void set_d_max(Vector const& d_max) { qp_.set_d_max(i_, d_max); }

			template <typename Vector>
			void set_x_min(Vector const& x_min) { qp_.set_x_min(i_, x_min - get_x()); }

			template <typename Vector>
			void set_x_max(Vector const& x_max) { qp_.set_x_max(i_, x_max - get_x()); }

			template <typename Vector>
			void set_u_min(Vector const& u_min) { qp_.set_u_min(i_, u_min - get_u()); }

			template <typename Vector>
			void set_u_max(Vector const& u_max) { qp_.set_u_max(i_, u_max - get_u()); }

		private:
			std::size_t const i_;
			QP& qp_;
			WorkingPoint workingPoint_;
			double const levenbergMarquardt_;
		};

		/**
		 * \brief Interface for setting MPC terminal stage data.
		 * OCP::InitTerminalStage() and OCP::UpdateTerminalStage() use this interface to set Hessian of the terminal stage cost,
		 * gradient of the terminal stage cost, linearized terminal stage constraints, and terminal state bounds.
		 *
		 * Hessian and gradient of Mayer term:
		 * H = [Q]
		 * g = [q]
		 *
		 */
		class TerminalStage
		{
		public:
			TerminalStage(WorkingPoint& wp, QP& qp, double lm)
			:	workingPoint_(wp), qp_(qp), levenbergMarquardt_(lm) {
				assert(wp.nT() == qp.nT());
			}

			TerminalStage(Stage const&) = delete;

			/**
			 * \brief Returns stage index.
			 */
			std::size_t get_index() const { return workingPoint_.nT(); }

			/**
			 * \brief Returns stage state.
			 */
			decltype(auto) get_x() const { return workingPoint_.get_x(get_index()); }

			template <typename Matrix>
			void set_Q(Matrix const &Q)	{
				qp_.set_Q(get_index(), Q + levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Vector>
			void set_q(Vector const &q) { qp_.set_q(get_index(), q); }

			template <typename Matrix>
			void set_C(Matrix const &C) { qp_.set_C_end(C); }

			template <typename Vector>
			void set_d_end_min(Vector const& d_min) { qp_.set_d_end_min(d_min); }

			template <typename Vector>
			void set_d_end_max(Vector const& d_max) { qp_.set_d_end_max(d_max); }

			template <typename Vector>
			void set_x_min(Vector const& x_min) { qp_.set_x_min(get_index(), x_min - get_x()); }

			template <typename Vector>
			void set_x_max(Vector const& x_max) { qp_.set_x_max(get_index(), x_max - get_x()); }

		private:
			QP& qp_;
			WorkingPoint workingPoint_;
			double const levenbergMarquardt_;
		};

		// Uses ocp_ to initialize _QP based on current working point workingPoint_.
		void InitQP()
		{
			for (unsigned i = 0; i < nT(); ++i)
			{
				Stage stage(i, workingPoint_, _QP, _levenbergMarquardt);
				_ocp.InitStage(stage);
			}

			TerminalStage terminal_stage(workingPoint_, _QP, _levenbergMarquardt);
			_ocp.InitTerminalStage(terminal_stage);
		}

		// Uses ocp_ to update _QP based on current working point workingPoint_.
		void UpdateQP()
		{
			for (unsigned i = 0; i < nT(); ++i)
			{
				Stage stage(i, workingPoint_, _QP, _levenbergMarquardt);
				_ocp.UpdateStage(stage);
			}

			TerminalStage terminal_stage(workingPoint_, _QP, _levenbergMarquardt);
			_ocp.UpdateTerminalStage(terminal_stage);
		}

		// Private data members.
		OCP const& _ocp;
		QP _QP;
		Solution solution_;
		QPSolver& solver_;

		double _levenbergMarquardt = 0.;

		/// \brief The working point (linearization point).
		WorkingPoint workingPoint_;

		// Preparation() sets this flag to true, Feedback() resets it to false.
		bool _prepared;
	};
}
