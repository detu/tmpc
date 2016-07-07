#pragma once

#include "Trajectory.hpp"
#include "OptimalControlProblem.hpp"
#include "matrix.hpp"

#include <stdexcept>
//#include <type_traits>

namespace tmpc
{
	template <typename OCP_, typename Integrator_, typename WorkingPoint, typename QP_>
	void realtime_iteration_preparation(OCP_ const& ocp, Integrator_ const& integrator,
			WorkingPoint const& working_point, QP_& qp, double levenberg_marquardt = 0.)
	{
		typedef typename OCP_::StateVector StateVector;
		typedef typename OCP_::InputVector InputVector;
		typedef typename OCP_::StateInputVector StateInputVector;
		typedef typename OCP_::ODEJacobianMatrix ODEJacobianMatrix;
		typedef typename OCP_::LagrangeHessianMatrix LagrangeHessianMatrix;
		typedef typename OCP_::MayerHessianMatrix MayerHessianMatrix;

		auto const N = ocp.getNumberOfIntervals();

		for (unsigned i = 0; i < N; ++i)
		{
			// Hessians and gradients of Lagrange terms.
			//
			LagrangeHessianMatrix H_i;
			StateInputVector g_i;

			StateInputVector z_i;
			z_i << working_point.get_x(i), working_point.get_u(i);

			ocp.LagrangeTerm(i, z_i, g_i, H_i);

			// Adding Levenberg-Marquardt term to make H positive-definite.
			set_H(qp, i, H_i + levenberg_marquardt * LagrangeHessianMatrix::Identity());
			set_g(qp, i, g_i);

			// Bound constraints.
			StateInputVector z_min, z_max;
			z_min << ocp.getStateMin(), ocp.getInputMin();
			z_max << ocp.getStateMax(), ocp.getInputMax();

			// C = [ssA, ssB];
			// x_{k+1} = C * z_k + c_k
			StateVector x_plus;
			ODEJacobianMatrix J;
			integrator.Integrate(i * integrator.timeStep(), z_i, x_plus, J);
			set_AB(qp, i, J);

			// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
			// c = f(z_k) - x_{k+1}
			qp.set_b(i, x_plus - working_point.get_x(i + 1));

			typename OCP_::ConstraintJacobianMatrix D;
			typename OCP_::ConstraintVector d_min, d_max;
			ocp.PathConstraints(i, z_i, D, d_min, d_max);
			set_CD(qp, i, D);
			qp.set_d_min(i, d_min);
			qp.set_d_max(i, d_max);

			// z_min stores _Nt vectors of size _Nz and 1 vector of size _Nx
			set_xu_min(qp, i, z_min - z_i);

			// z_max stores _Nt vectors of size _Nz and 1 vector of size _Nx
			set_xu_max(qp, i, z_max - z_i);
		}

		auto const xN = working_point.get_x(N);

		// End state constraints.
		typename OCP_::TerminalConstraintJacobianMatrix D;
		typename OCP_::TerminalConstraintVector d_min, d_max;
		ocp.TerminalConstraints(xN, D, d_min, d_max);
		qp.set_C_end(D);
		qp.set_d_end_min(d_min);
		qp.set_d_end_max(d_max);

		qp.set_x_min(N, ocp.getTerminalStateMin() - xN);
		qp.set_x_max(N, ocp.getTerminalStateMax() - xN);

		// Hessian and gradient of Mayer term.
		typename OCP_::MayerHessianMatrix H_T;
		typename OCP_::StateVector g_T;
		ocp.MayerTerm(xN, g_T, H_T);

		// Adding Levenberg-Marquardt term to make H positive-definite.
		qp.set_Q(N, H_T + levenberg_marquardt * MayerHessianMatrix::Identity());
		qp.set_q(N, g_T);

	}

	template <typename StateVector, typename WorkingPoint, typename QP, typename QPSolver, typename Solution>
	decltype(auto) realtime_iteration_feedback(StateVector& x0, WorkingPoint& working_point, QP& qp, QPSolver& solver, Solution& solution)
	{
		/** embed current initial value */
		auto const w0 = (x0 - working_point.get_x(0)).eval();
		qp.set_x_min(0, w0);
		qp.set_x_max(0, w0);

		/** solve QP */
		solver.Solve(qp, solution);

		// Add QP step to the working point.
		working_point += solution;

		// Return the calculated control input.
		return working_point.get_u(0);
	}

	template<class _Problem, typename Integrator_, class QPSolver_>
	class RealtimeIteration
	{
	public:
		typedef _Problem Problem;
		typedef QPSolver_ QPSolver;
		typedef Integrator_ Integrator;

		typedef typename QPSolver::Solution Solution;
		typedef Trajectory<Problem::NX, Problem::NU> WorkingPoint;

		typedef typename Problem::StateVector StateVector;
		typedef typename Problem::InputVector InputVector;
		typedef typename Problem::StateInputVector StateInputVector;
		typedef typename Problem::ODEJacobianMatrix ODEJacobianMatrix;
		typedef typename Problem::LagrangeHessianMatrix LagrangeHessianMatrix;
		typedef typename Problem::MayerHessianMatrix MayerHessianMatrix;

		typedef std::function<void (typename QPSolver::Problem const&)> QPCallback;

		RealtimeIteration(Problem const& ocp, Integrator const& integrator, QPSolver& solver, WorkingPoint const& working_point)
		:	_ocp(ocp)
		,	_QP(working_point.nT())
		,	_workingPoint(working_point)
		,	_solution(working_point.nT())
		,	_Solver(solver)
		,	_levenbergMarquardt(0.0)
		,	_integrator(integrator)
		,	_prepared(true)
		{
			UpdateQP();
		}

		RealtimeIteration(RealtimeIteration const&) = delete;

		// Feed current state x0, get back control input u.
		InputVector Feedback(const StateVector& x0)
		{
			if (!_prepared)
				throw std::logic_error("ModelPredictiveController::Feedback(): controller is not prepared.");

			auto const u = realtime_iteration_feedback(x0, _workingPoint, _QP, _Solver, _solution);

			_prepared = false;

			// Return the calculated control input.
			return u;
		}

		void Preparation()
		{
			if (_prepared)
				throw std::logic_error("ModelPredictiveController::Preparation(): controller is already prepared.");

			shift(_workingPoint);

			// Calculate new QP.
			UpdateQP();

			_prepared = true;
		}

		WorkingPoint const& getWorkingPoint() const { return _workingPoint; }

		double getLevenbergMarquardt() const { return _levenbergMarquardt; }
		void setLevenbergMarquardt(double val) { _levenbergMarquardt = val; }

		unsigned nT() const { return _ocp.getNumberOfIntervals(); }
		unsigned nU() const	noexcept { return _Nu; }
		unsigned nX() const	noexcept { return _Nx; }
		unsigned nZ() const noexcept { return _Nz; }

	private:

		// Initializes _G, _g, _y, _C, _c, _zMin, _zMax based on current working point _w.
		void UpdateQP()
		{
			realtime_iteration_preparation(_ocp, _integrator, _workingPoint, _QP, _levenbergMarquardt);
		}

		// Private data members.
		Problem const& _ocp;
		Integrator const& _integrator;

		static const unsigned _Nu = Problem::NU;
		static const unsigned _Nx = Problem::NX;
		static const unsigned _Nz = Problem::NW;
		static const unsigned _Nd = Problem::NC;
		static const unsigned _NdT = Problem::NCT;
		
		typename QPSolver::Problem _QP;
		typename QPSolver::Solution _solution;
		QPSolver& _Solver;

		double _levenbergMarquardt;

		// Working point (linearization point).
		WorkingPoint _workingPoint;

		// Preparation() sets this flag to true, Feedback() resets it to false.
		bool _prepared;
	};
}
