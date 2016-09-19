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
	template <typename Problem, typename Integrator_,
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
		typedef Integrator_ Integrator;
		typedef Trajectory<NX, NU, NW, NY> WorkingPoint;

		/**
		 * \brief Constructor.
		 * TODO: Change aggregation to containment and init using rvalue references?
		 */
		MovingHorizonEstimator(Problem const& problem, Integrator const& integrator,
				QPSolver& solver, WorkingPoint const& working_point)
		:	problem_(problem)
		,	_QP(working_point.nT())
		,	workingPoint_(working_point)
		,	_solution(working_point.nT())
		,	solver_(solver)
		,	_integrator(integrator)
		,	_prepared(true)
		{
			UpdateQP();
		}

		MovingHorizonEstimator(MovingHorizonEstimator const&) = delete;

		/**
		 * \brief Number of intervals.
		 */
		std::size_t nT() const { return workingPoint_.nT(); }

		/**
		 * \brief Feed input u[N-1] and measurement y[N-1], get back current state estimate x[N].
		 * \param [in] u -- input vector at previous time instant
		 * \param [out] y -- measurement vector at previous time instant
		 * \return State estimate at current time instant
		 */
		template <typename InputVector, typename OutputVector>
		decltype(auto) Feedback(InputVector const& u, OutputVector const& y)
		{
			if (!_prepared)
				throw std::logic_error("MovingHorizonEstimator::Feedback(): estimator is not prepared.");

			// Update the last stage of QP.
			workingPoint_.set_u(nT() - 1, u);
			workingPoint_.set_y(nT() - 1, y);
			UpdateQP(nT() - 1);

			/** solve QP */
			solver_.Solve(_QP, _solution);

			// Add QP step to the working point.
			{
				for (std::size_t i = 0; i < nT() + 1; ++i)
					workingPoint_.set_x(i, workingPoint_.get_x(i) + _solution.get_x(i));

				for (std::size_t i = 0; i < nT(); ++i)
					workingPoint_.set_w(i, workingPoint_.get_w(i) + _solution.get_u(i));
			}

			_prepared = false;

			// Return the estimated state.
			return workingPoint_.get_x(nT());
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

		// Initializes _G, _g, _y, _C, _c, _zMin, _zMax based on current working point _w.
		void UpdateQP()
		{
			auto const N = nT();

			for (unsigned i = 0; i < N; ++i)
			{
				// Hessian and gradient of stage cost.
				// H = [Q, S
				//      S, R]
				// g = [q
				//      r]

				Eigen::Matrix<double, NX, NX> Q;
				Eigen::Matrix<double, NU, NU> R;
				Eigen::Matrix<double, NX, NU> S;
				Eigen::Matrix<double, NX,  1> q;
				Eigen::Matrix<double, NU,  1> r;

				// TODO: change the interface of QP so that it can one can write
				// _ocp.LagrangeTerm(... , ... , ... , qp.Q(i), qp.R(i), qp.S(i), qp.q(i), qp.r(i));
				// i.e. Q(i), R(i), etc. are references of some proxy objects?
				// HINT: Use rvalue references?
				problem_.StageCost(i, workingPoint_.get_x(i), workingPoint_.get_u(i),
						              workingPoint_.get_w(i), workingPoint_.get_y(i), Q, R, S, q, r);

				_QP.set_Q(i, Q + _levenbergMarquardt * identity<decltype(Q)>());	// Adding Levenberg-Marquardt term to make H positive-definite.
				_QP.set_R(i, R + _levenbergMarquardt * identity<decltype(R)>());	// Adding Levenberg-Marquardt term to make H positive-definite.
				_QP.set_S(i, S);
				_QP.set_q(i, q);
				_QP.set_r(i, r);

				// C = [ssA, ssB];
				// x_{k+1} = C * z_k + c_k
				Eigen::Matrix<double, NX, 1> x_plus;
				Eigen::Matrix<double, NX, NX> A;
				Eigen::Matrix<double, NX, NU> B;
				integrate(_integrator, problem_.getODE(workingPoint_.get_u(i)),
						i * _integrator.timeStep(), workingPoint_.get_x(i), workingPoint_.get_w(i), x_plus, A, B);
				_QP.set_A(i, A);
				_QP.set_B(i, B);

				// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
				// c = f(z_k) - x_{k+1}
				_QP.set_b(i, x_plus - workingPoint_.get_x(i + 1));

				Eigen::Matrix<double, NC, NX + NU> D;
				Eigen::Matrix<double, NC, 1> d_min, d_max;
				
				{
					Eigen::Matrix<double, NX + NU, 1> z_i;
					z_i << workingPoint_.get_x(i), workingPoint_.get_u(i);
					problem_.PathConstraints(i, z_i, D, d_min, d_max);
				}

				set_CD(_QP, i, D);
				_QP.set_d_min(i, d_min);
				_QP.set_d_max(i, d_max);
				
				// Bound constraints.
				_QP.set_x_min(i, problem_.getStateMin() - workingPoint_.get_x(i));	_QP.set_u_min(i, problem_.getInputMin() - workingPoint_.get_u(i));
				_QP.set_x_max(i, problem_.getStateMax() - workingPoint_.get_x(i));	_QP.set_u_max(i, problem_.getInputMax() - workingPoint_.get_u(i));
			}

			auto const xN = workingPoint_.get_x(N);

			// Hessian and gradient of arrival cost.
			Eigen::Matrix<double, NX, NX> H_arr;
			Eigen::Matrix<double, NX, 1> g_arr;
			problem_.ArrivalCost(workingPoint_.get_x(0), g_arr, H_arr);

			// Arrival cost is added to the cost of stage 0.
			_QP.set_Q(0, _QP.get_Q(0) + H_arr);
			_QP.set_q(0, _QP.get_q(0) + g_arr);

			// The total cost does not depend on x[N],
			// but we add Levenberg-Marquardt term to H[N] make H positive-definite.
			_QP.set_Q(N, _levenbergMarquardt * identity<Eigen::Matrix<double, NX, NX>>());
			_QP.set_q(N, zero<Eigen::Matrix<double, NX, 1>>());

			// Terminal state constraints are empty.
			assert(_QP.nDT() == 0);

			// Bounds for terminal state the same as for other states.
			_QP.set_x_min(N, problem_.getStateMin() - xN);
			_QP.set_x_max(N, problem_.getStateMax() - xN);
		}

		// Private data members.

		Problem const& problem_;
		Integrator const& _integrator;
		
		typename QPSolver::Problem _QP;
		typename QPSolver::Solution _solution;
		QPSolver& solver_;

		double _levenbergMarquardt = 0.;

		// Working point (linearization point).
		WorkingPoint workingPoint_;

		// Preparation() sets this flag to true, Feedback() resets it to false.
		bool _prepared;
	};
}
