#pragma once

#include "Trajectory.hpp"
#include "matrix.hpp"
#include "../qp/qp.hpp"

#include <stdexcept>
//#include <type_traits>

namespace tmpc
{
	/**
	 * \brief Implements moving horizon estimation algorithm.
	 * \tparam <Problem> A class describing Trajectory Estimation Problem.
	 */
	template <typename Problem,
		template <unsigned, unsigned, unsigned, unsigned> class QPSolver_>
	class MovingHorizonEstimator
	{
		static auto const NX = Problem::NX;
		static auto const NU = Problem::NU;
		static auto const NW = Problem::NW;
		static auto const NY = Problem::NY;
		static auto const NC = Problem::NC;

	public:
		typedef QPSolver_<NX, NW, NC, 0> QPSolver;
		typedef Trajectory<NX, NU, NW, NY> WorkingPoint;

		/**
		 * \brief Constructor.
		 * TODO: Change aggregation to containment and init using rvalue references?
		 */
		MovingHorizonEstimator(Problem const& problem,
				QPSolver& solver, WorkingPoint const& working_point)
		:	problem_(problem)
		,	solution_(working_point.nT())
		,	solver_(solver)
		,	work_(working_point)
		,	_prepared(true)
		{
			UpdateQP();
		}

		MovingHorizonEstimator(MovingHorizonEstimator const&) = delete;

		/**
		 * \brief Number of intervals.
		 */
		std::size_t nT() const { return work_.workingPoint_.nT(); }

		/**
		 * \brief Feed input u[k] and measurement y[k], get back current state estimate x[k+1].
		 * \param [in] u -- input vector at time interval k
		 * \param [out] y -- measurement vector at time instant k
		 * \return State estimate at time instant k+1
		 */
		template <typename InputVector, typename OutputVector>
		decltype(auto) Feedback(InputVector const& u, OutputVector const& y)
		{
			if (!_prepared)
				throw std::logic_error("MovingHorizonEstimator::Feedback(): estimator is not prepared.");

			// Update the last stage of the QP.
			work_.workingPoint_.set_u(nT() - 1, u);
			work_.workingPoint_.set_y(nT() - 1, y);
			UpdateQP(nT() - 1);

			/** solve QP */
			solver_.Solve(work_.qp_, solution_);

			// Add QP step to the working point.
			{
				for (std::size_t i = 0; i < nT() + 1; ++i)
					work_.workingPoint_.set_x(i, work_.workingPoint_.get_x(i) + solution_.get_x(i));

				for (std::size_t i = 0; i < nT(); ++i)
					work_.workingPoint_.set_w(i, work_.workingPoint_.get_w(i) + solution_.get_u(i));
			}

			_prepared = false;

			// Return the estimated state.
			return work_.workingPoint_.get_x(nT());
		}

		void Preparation()
		{
			if (_prepared)
				throw std::logic_error("MovingHorizonEstimator::Preparation(): estimator is already prepared.");

			shift(work_.workingPoint_);

			// Calculate new QP.
			UpdateQP();

			_prepared = true;
		}

		WorkingPoint const& getWorkingPoint() const { return work_.workingPoint_; }

		double getLevenbergMarquardt() const { return work_.levenbergMarquardt_; }
		void setLevenbergMarquardt(double val) { work_.levenbergMarquardt_ = val; }

		/*
		unsigned nT() const { return _ocp.getNumberOfIntervals(); }
		unsigned nU() const	noexcept { return _Nu; }
		unsigned nX() const	noexcept { return _Nx; }
		unsigned nZ() const noexcept { return _Nz; }
		*/

	private:
		typedef typename QPSolver::Solution Solution;
		typedef typename QPSolver::Problem QP;

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
			 * NOTICE that "u" in QP corresponds to the disturbance "w".
			 */
			~Stage()
			{
				work_.qp_.set_x_min(i_, work_.lowerBound_.get_x(i_) - get_x());
				work_.qp_.set_x_max(i_, work_.upperBound_.get_x(i_) - get_x());
				work_.qp_.set_u_min(i_, work_.lowerBound_.get_w(i_) - get_w());
				work_.qp_.set_u_max(i_, work_.upperBound_.get_w(i_) - get_w());
			}

			template <typename Matrix>
			void set_Q(Matrix const &Q)	{
				work_.qp_.set_Q(i_, Q + work_.levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_R(Matrix const &R) {
				work_.qp_.set_R(i_, R + work_.levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Matrix>
			void set_S(Matrix const &S) { work_.qp_.set_S(i_, S); }

			template <typename Vector>
			void set_q(Vector const &q) { work_.qp_.set_q(i_, q); }

			template <typename Vector>
			void set_r(Vector const &r) { work_.qp_.set_r(i_, r); }

			template <typename Matrix>
			void set_A(Matrix const &A) { work_.qp_.set_A(i_, A); }

			template <typename Matrix>
			void set_B(Matrix const &B) { work_.qp_.set_B(i_, B); }

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
			void set_D(Matrix const &D) { work_.qp_.set_D(i_, D); }

			template <typename Vector>
			void set_d_min(Vector const& d_min) { work_.qp_.set_d_min(i_, d_min); }

			template <typename Vector>
			void set_d_max(Vector const& d_max) { work_.qp_.set_d_max(i_, d_max); }

			template <typename Vector>
			void set_x_min(Vector const& x_min) { work_.lowerBound_.set_x(i_, x_min); }

			template <typename Vector>
			void set_x_max(Vector const& x_max) { work_.upperBound_.set_x(i_, x_max); }

			template <typename Vector>
			void set_u_min(Vector const& u_min) { work_.lowerBound_.set_u(i_, u_min); }

			template <typename Vector>
			void set_u_max(Vector const& u_max) { work_.upperBound_.set_u(i_, u_max); }

		private:
			std::size_t const i_;
			Work& work_;

			decltype(auto) get_x() const { return work_.workingPoint_.get_x(i_); }
			decltype(auto) get_w() const { return work_.workingPoint_.get_w(i_); }
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

			template <typename Matrix>
			void set_Q(Matrix const &Q)	{
				work_.qp_.set_Q(get_index(), Q + work_.levenbergMarquardt_ * identity<Matrix>());	// Adding Levenberg-Marquardt term to make H positive-definite.
			}

			template <typename Vector>
			void set_q(Vector const &q) { work_.qp_.set_q(get_index(), q); }

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

			decltype(auto) get_x() const { return work_.workingPoint_.get_x(get_index()); }
			decltype(auto) get_index() const { return work_.workingPoint_.nT(); }
		};

		// Uses ocp_ to initialize _QP based on current working point workingPoint_.
		void InitQP()
		{
			for (unsigned i = 0; i < nT(); ++i)
			{
				Stage stage(i, work_);
				problem_.InitStage(i, stage);
			}

			TerminalStage terminal_stage(work_);
			problem_.InitTerminalStage(terminal_stage);
		}

		/// \brief Internal work data structure.
		struct Work
		{
			double levenbergMarquardt_ = 0.;

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
					constant<typename WorkingPoint::StateVector>(-std::numeric_limits<double>::infinity()),
					constant<typename WorkingPoint::InputVector>(-std::numeric_limits<double>::infinity()))
			,	upperBound_(working_point.nT(),
					constant<typename WorkingPoint::StateVector>( std::numeric_limits<double>::infinity()),
					constant<typename WorkingPoint::InputVector>( std::numeric_limits<double>::infinity()))
			,	qp_(working_point.nT())
			{
			}
		};

		// Uses ocp_ to update _QP based on current working point workingPoint_.
		void UpdateQP()
		{
			for (unsigned i = 0; i < nT(); ++i)
				UpdateQP(i);

			TerminalStage terminal_stage(work_);
			problem_.UpdateTerminalStage(work_.workingPoint_.get_x(nT()), terminal_stage);
		}

		void UpdateQP(unsigned i)
		{
			assert(i < nT());

			Stage stage(i, work_);
			problem_.UpdateStage(i, work_.workingPoint_.get_x(i), work_.workingPoint_.get_u(i),
					work_.workingPoint_.get_w(i), work_.workingPoint_.get_y(i), stage);
		}

		// Private data members.
		Problem const& problem_;
		Solution solution_;
		QPSolver& solver_;
		Work work_;

		// Preparation() sets this flag to true, Feedback() resets it to false.
		bool _prepared;
	};
}
