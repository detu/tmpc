#pragma once

#include "../qp/diagnostics.hpp"
#include "Trajectory.hpp"
#include "../qp/qp.hpp"

#include <stdexcept>
#include <limits>
#include <iterator>
//#include <type_traits>

namespace tmpc
{
	/**
	 * \brief Implements an MPC controller with realtime iteration scheme.
	 *
	 * \tparam <K> Class implementing the Kernel concept
	 * \tparam <D> Class defining problem dimensions
	 */
	template <typename K, typename D, typename OCP, typename QPSolver_>
	class RealtimeIteration
	{
		static auto constexpr NX = D::NX;
		static auto constexpr NU = D::NU;
		static auto constexpr NC = D::NC;
		static auto constexpr NCT = D::NCT;

	public:
		typedef QPSolver_ QPSolver;
		typedef typename QPSolver::Problem QP;
		typedef typename QPSolver::Solution Solution;

		// TODO: parameterize Trajectory with K
		typedef Trajectory<NX, NU> WorkingPoint;

		RealtimeIteration(OCP const& ocp, QPSolver& solver, WorkingPoint const& working_point)
		:	_ocp(ocp)
		,	solution_(working_point.nT())
		,	solver_(solver)
		,	work_(working_point)
		,	_prepared(false)
		{
			InitQP();
		}

		RealtimeIteration(RealtimeIteration const&) = delete;

		/**
		 * \brief Number of intervals.
		 */
		std::size_t nT() const { return work_.workingPoint_.nT(); }

		// Feed current state x0, get back control input u.
		template <typename StateVector>
		decltype(auto) Feedback(const StateVector& x0)
		{
			if (!_prepared)
				throw std::logic_error("RealtimeIteration::Feedback(): RTI is not prepared.");

			/** embed current initial value */
			auto const w0 = K::eval(x0 - work_.workingPoint_.get_x(0));
			work_.qp_.set_x_min(0, w0);
			work_.qp_.set_x_max(0, w0);

			/** solve QP */
			solver_.Solve(work_.qp_, solution_);

			// Add QP step to the working point.
			{
				for (std::size_t i = 0; i < nT() + 1; ++i)
					work_.workingPoint_.set_x(i, work_.workingPoint_.get_x(i) + solution_.get_x(i));

				for (std::size_t i = 0; i < nT(); ++i)
					work_.workingPoint_.set_u(i, work_.workingPoint_.get_u(i) + solution_.get_u(i));
			}

			_prepared = false;

			// Return the calculated control input.
			return work_.workingPoint_.get_u(0);
		}

		/// \brief Number of iterations performed by the QP solver during the last Feedback phase.
		unsigned getNumIter() const { return solution_.getNumIter(); }

		void Preparation()
		{
			if (_prepared)
				throw std::logic_error("RealtimeIteration::Preparation(): RTI is already prepared.");

			// MK: Working point shifting is disabled. The main reason is that it is possible to have
			// time step in MPC different (typicaly bigger than) from the controller sampling time.
			// Perhaps, shifting policy can be defined externally.

			// shift(work_.workingPoint_);

			// Calculate new QP.
			UpdateQP();

			_prepared = true;
		}

		WorkingPoint const& getWorkingPoint() const { return work_.workingPoint_; }

		// TODO: change LevenbergMarquardt to be an argument of Preparation()?
		typename K::Scalar getLevenbergMarquardt() const { return work_.levenbergMarquardt_; }

		void setLevenbergMarquardt(typename K::Scalar val)
		{
			if (_prepared)
				throw std::logic_error("RealtimeIteration::setLevenbergMarquardt(): RTI must not be prepared when setting LevenbergMarquardt.");

			if (val < 0.)
				throw std::invalid_argument("RealtimeIteration::setLevenbergMarquardt(): LevenbergMarquardt must be non-negative.");

			if (work_.levenbergMarquardt_ != val)
			{
				work_.levenbergMarquardt_ = val;
				InitQP();	// Re-initialize QPs because of the changed levenbergMarquardt_ value.
			}
		}

		QP const& getQP() const
		{
			return work_.qp_;
		}

		Solution const& getSolution() const
		{
			return solution_;
		}

		/*
		unsigned nT() const { return _ocp.getNumberOfIntervals(); }
		unsigned nU() const	noexcept { return _Nu; }
		unsigned nX() const	noexcept { return _Nx; }
		unsigned nZ() const noexcept { return _Nz; }
		*/

	private:

		struct Work;

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
			Stage(std::size_t i, Work& work)
			:	i_(i), work_(work) {
			}

			Stage(Stage const&) = delete;

			/**
			 * \brief On destruction, updates bounds of the QP.
			 */
			~Stage()
			{
				work_.qp_.set_x_min(i_, work_.lowerBound_.get_x(i_) - get_x());
				work_.qp_.set_x_max(i_, work_.upperBound_.get_x(i_) - get_x());
				work_.qp_.set_u_min(i_, work_.lowerBound_.get_u(i_) - get_u());
				work_.qp_.set_u_max(i_, work_.upperBound_.get_u(i_) - get_u());
			}

			/**
			 * \brief Returns stage index.
			 */
			std::size_t get_index() const { return i_; }

			/**
			 * \brief Returns stage state.
			 */
			decltype(auto) get_x() const { return work_.workingPoint_.get_x(i_); }

			/**
			 * \brief Returns stage input.
			 */
			decltype(auto) get_u() const { return work_.workingPoint_.get_u(i_); }

			template <typename Matrix>
			void set_Q(Matrix const &Q)	{
				work_.qp_.set_Q(i_, Q + work_.levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_R(Matrix const &R) {
				work_.qp_.set_R(i_, R + work_.levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_R(typename K::size_t i, typename K::size_t j, Matrix const &R)
			{
				auto constexpr M = K::template rows<Matrix>();
				auto constexpr N = K::template cols<Matrix>();

				auto const I = K::template identity<D::NU, D::NU>();
				work_.qp_.set_R(i_, i, j, R + work_.levenbergMarquardt_ * K::template block<M, N>(I, i, j));
										   // ^ Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_S(Matrix const &S) { work_.qp_.set_S(i_, S); }

			template <typename Matrix>
			void set_S(typename K::size_t i, typename K::size_t j, Matrix const &S)
			{
				work_.qp_.set_S(i_, i, j, S);
			}

			template <typename Vector>
			void set_q(Vector const &q) { work_.qp_.set_q(i_, q); }

			template <typename Vector>
			void set_r(Vector const &r) { work_.qp_.set_r(i_, r); }

			template <typename Vector>
			void set_r(typename K::size_t i, Vector const &r)
			{
				work_.qp_.set_r(i_, i, r);
			}

			template <typename Matrix>
			void set_A(Matrix const &A) { work_.qp_.set_A(i_, A); }

			template <typename Matrix>
			void set_B(Matrix const &B) { work_.qp_.set_B(i_, B); }

			template <typename Matrix>
			void set_B(typename K::size_t i, typename K::size_t j, Matrix const &B)
			{
				work_.qp_.set_B(i_, i, j, B);
			}

			template <typename Vector>
			void set_x_next(Vector const& x_next)
			{
				// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
				// c = f(z_k) - x_{k+1}
				work_.qp_.set_b(i_, x_next - work_.workingPoint_.get_x(i_ + 1));
			}

			template <typename Matrix>
			void set_C(Matrix const &C) { work_.qp_.set_C(i_, C); }

			template <typename Matrix>
			void set_D(Matrix const &val) { work_.qp_.set_D(i_, val); }

			template <typename Vector>
			void set_d_min(Vector const& d_min) { work_.qp_.set_d_min(i_, d_min); }

			template <typename Vector>
			void set_d_max(Vector const& d_max) { work_.qp_.set_d_max(i_, d_max); }

			template <typename Vector>
			void set_x_min(Vector const& x_min) { work_.lowerBound_.set_x(i_, x_min); }

			template <typename Vector>
			void set_x_max(Vector const& x_max) { work_.upperBound_.set_x(i_, x_max); }

			template <typename Vector>
			void set_u_min(Vector const& u_min)
			{
				work_.lowerBound_.set_u(i_, u_min);
			}

			template <typename Vector>
			void set_u_min(typename K::size_t i, Vector const& u_min)
			{
				work_.lowerBound_.set_u(i_, i, u_min);
			}

			template <typename Vector>
			void set_u_max(Vector const& u_max)
			{
				work_.upperBound_.set_u(i_, u_max);
			}

			template <typename Vector>
			void set_u_max(typename K::size_t i, Vector const& u_max)
			{
				work_.upperBound_.set_u(i_, i, u_max);
			}

		private:
			std::size_t const i_;
			Work& work_;
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
			TerminalStage(Work& work)
			:	work_(work)
			{
			}

			TerminalStage(Stage const&) = delete;

			/**
			 * \brief On destruction, updates bounds of the QP.
			 */
			~TerminalStage()
			{
				work_.qp_.set_x_min(get_index(), work_.lowerBound_.get_x(get_index()) - get_x());
				work_.qp_.set_x_max(get_index(), work_.upperBound_.get_x(get_index()) - get_x());
			}

			/**
			 * \brief Returns stage index.
			 */
			std::size_t get_index() const { return work_.workingPoint_.nT(); }

			/**
			 * \brief Returns stage state.
			 */
			decltype(auto) get_x() const { return work_.workingPoint_.get_x(get_index()); }

			// TODO: Isn't it more flexible and clear to add Lev-Mar in the OCP?
			template <typename Matrix>
			void set_Q(Matrix const &Q)	{
				work_.qp_.set_Q(get_index(), Q + work_.levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			decltype(auto) get_Q() const
			{
				return work_.qp_.get_Q(get_index());
			}

			template <typename Vector>
			void set_q(Vector const &q)
			{
				work_.qp_.set_q(get_index(), q);
			}

			decltype(auto) get_q() const
			{
				return work_.qp_.get_q(get_index());
			}

			template <typename Matrix>
			void set_C(Matrix const &C) { work_.qp_.set_C_end(C); }

			template <typename Vector>
			void set_d_end_min(Vector const& d_min) { work_.qp_.set_d_end_min(d_min); }

			template <typename Vector>
			void set_d_end_max(Vector const& d_max) { work_.qp_.set_d_end_max(d_max); }

			template <typename Vector>
			void set_x_min(Vector const& x_min) { work_.lowerBound_.set_x(get_index(), x_min); }

			template <typename Vector>
			void set_x_max(Vector const& x_max) { work_.upperBound_.set_x(get_index(), x_max); }

		private:
			Work& work_;
		};

		/**
		 * \brief Set all data of the QP to NaN.
		 */
		void SetQpNaN()
		{
			auto constexpr nan = std::numeric_limits<typename K::Scalar>::signaling_NaN();
			auto& qp = work_.qp_;

			for (unsigned i = 0; i < nT(); ++i)
			{
				qp.set_Q(i, K::template constant<D::NX, D::NX>(nan));
				qp.set_R(i, K::template constant<D::NU, D::NU>(nan));
				qp.set_S(i, K::template constant<D::NX, D::NU>(nan));
				qp.set_q(i, K::template constant<D::NX>(nan));
				qp.set_r(i, K::template constant<D::NU>(nan));

				qp.set_A(i, K::template constant<D::NX, D::NX>(nan));
				qp.set_B(i, K::template constant<D::NX, D::NU>(nan));

				qp.set_x_min(i, K::template constant<D::NX>(nan));
				qp.set_x_max(i, K::template constant<D::NX>(nan));
				qp.set_u_min(i, K::template constant<D::NU>(nan));
				qp.set_u_max(i, K::template constant<D::NU>(nan));

				qp.set_C(i, K::template constant<D::NC, D::NX>(nan));
				qp.set_D(i, K::template constant<D::NC, D::NU>(nan));

				qp.set_d_min(i, K::template constant<D::NC>(nan));
				qp.set_d_max(i, K::template constant<D::NC>(nan));
			}

			set_Q_end(qp, K::template constant<D::NX, D::NX>(nan));
			set_q_end(qp, K::template constant<D::NX>(nan));

			set_x_end_min(qp, K::template constant<D::NX>(nan));
			set_x_end_max(qp, K::template constant<D::NX>(nan));

			qp.set_C_end(K::template constant<D::NCT, D::NX>(nan));
			qp.set_d_end_min(K::template constant<D::NCT>(nan));
			qp.set_d_end_max(K::template constant<D::NCT>(nan));
		}

		// Uses ocp_ to initialize _QP based on current working point workingPoint_.
		void InitQP()
		{
			SetQpNaN();

			for (unsigned i = 0; i < nT(); ++i)
			{
				Stage stage(i, work_);
				_ocp.InitStage(stage);
			}

			TerminalStage terminal_stage(work_);
			_ocp.InitTerminalStage(terminal_stage);
		}

		// Uses ocp_ to update _QP based on current working point workingPoint_.
		void UpdateQP()
		{
			for (unsigned i = 0; i < nT(); ++i)
			{
				Stage stage(i, work_);
				_ocp.UpdateStage(stage);
			}

			TerminalStage terminal_stage(work_);
			_ocp.UpdateTerminalStage(terminal_stage);

#ifndef NDEBUG
			// Error catching: check if we have any unexpected NaNs of infs in QP after update.

			std::vector<std::string> messages;
			qp::diagnose<K>(work_.qp_, std::back_inserter(messages));

			if (!messages.empty())
			{
				std::stringstream err_msg;
				err_msg << "RealtimeIteration encountered an ill-formed QP. Messages follow:" << std::endl;
				std::copy(messages.begin(), messages.end(), std::ostream_iterator<std::string>(err_msg, "\n"));

				throw std::logic_error(err_msg.str());
			}
#endif
		}

		/// \brief Internal work data structure.
		struct Work
		{
			typename K::Scalar levenbergMarquardt_ = 0.;

			/// \brief The working point (linearization point).
			WorkingPoint workingPoint_;

			/// \brief Lower bound of the working point.
			WorkingPoint lowerBound_;

			/// \brief Upper bound of the working point.
			WorkingPoint upperBound_;

			/// Multistage QP.
			QP qp_;

			Work(WorkingPoint const& working_point)
			:	workingPoint_(working_point)
			,	lowerBound_(working_point.nT(),
					K::template constant<D::NX>(-inf_),
					K::template constant<D::NU>(-inf_))
			,	upperBound_(working_point.nT(),
					K::template constant<D::NX>( inf_),
					K::template constant<D::NU>( inf_))
			,	qp_(working_point.nT())
			{
			}
		};

		// Private data members.
		OCP const& _ocp;	// <-- TODO: change to value semantics
		Solution solution_;
		QPSolver& solver_;	// <-- TODO: change to value semantics
		Work work_;

		// Preparation() sets this flag to true, Feedback() resets it to false.
		bool _prepared;

		static auto constexpr inf_ = std::numeric_limits<typename K::Scalar>::infinity();
	};
}
