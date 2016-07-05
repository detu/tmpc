#pragma once

#include "Trajectory.hpp"

#include <stdexcept>

namespace tmpc
{
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

			/** embed current initial value */
			auto const w0 = (x0 - _workingPoint.get_x(0)).eval();
			_QP.set_x_min(0, w0);
			_QP.set_x_max(0, w0);

			/** solve QP */
			_Solver.Solve(_QP, _solution);

			_prepared = false;

			// Return the calculated control input.
			return _workingPoint.get_u(0) + _solution.get_u(0);
		}

		void Preparation()
		{
			if (_prepared)
				throw std::logic_error("ModelPredictiveController::Preparation(): controller is already prepared.");

			// Add QP step to the working point.
			_workingPoint += _solution;

			/** prepare QP for next solution */
			//qpDUNES_shiftLambda(&_qpData);			/* shift multipliers */
			//qpDUNES_shiftIntervals(&_qpData);		/* shift intervals (particularly important when using qpOASES for underlying local QPs) */

			// Shift working point
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
			auto const N = _ocp.getNumberOfIntervals();
			for (unsigned i = 0; i < N; ++i)
			{
				// Hessians and gradients of Lagrange terms.
				//
				LagrangeHessianMatrix H_i;
				StateInputVector g_i;

				StateInputVector z_i;
				z_i << _workingPoint.get_x(i), _workingPoint.get_u(i);

				_ocp.LagrangeTerm(i, z_i, g_i, H_i);

				// Adding Levenberg-Marquardt term to make H positive-definite.
				set_H(_QP, i, H_i + _levenbergMarquardt * LagrangeHessianMatrix::Identity());
				set_g(_QP, i, g_i);

				// Bound constraints.
				StateInputVector z_min, z_max;
				z_min << _ocp.getStateMin(), _ocp.getInputMin();
				z_max << _ocp.getStateMax(), _ocp.getInputMax();

				// C = [ssA, ssB];
				// x_{k+1} = C * z_k + c_k
				typename Problem::StateVector x_plus;
				typename Problem::ODEJacobianMatrix J;
				_integrator.Integrate(i * _integrator.timeStep(), z_i, x_plus, J);
				set_AB(_QP, i, J);

				// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
				// c = f(z_k) - x_{k+1}
				_QP.set_b(i, x_plus - _workingPoint.get_x(i + 1));

				typename Problem::ConstraintJacobianMatrix D;
				typename Problem::ConstraintVector d_min, d_max;
				_ocp.PathConstraints(i, z_i, D, d_min, d_max);
				set_CD(_QP, i, D);
				_QP.set_d_min(i, d_min);
				_QP.set_d_max(i, d_max);

				// z_min stores _Nt vectors of size _Nz and 1 vector of size _Nx
				set_xu_min(_QP, i, z_min - z_i);

				// z_max stores _Nt vectors of size _Nz and 1 vector of size _Nx
				set_xu_max(_QP, i, z_max - z_i);
			}

			auto const xN = _workingPoint.get_x(N);

			// End state constraints.
			typename Problem::TerminalConstraintJacobianMatrix D;
			typename Problem::TerminalConstraintVector d_min, d_max;
			_ocp.TerminalConstraints(xN, D, d_min, d_max);
			_QP.set_C_end(D);
			_QP.set_d_end_min(d_min);
			_QP.set_d_end_max(d_max);

			_QP.set_x_min(N, _ocp.getTerminalStateMin() - xN);
			_QP.set_x_max(N, _ocp.getTerminalStateMax() - xN);

			// Hessian and gradient of Mayer term.
			typename Problem::MayerHessianMatrix H_T;
			typename Problem::StateVector g_T;
			_ocp.MayerTerm(xN, g_T, H_T);

			// Adding Levenberg-Marquardt term to make H positive-definite.
			_QP.set_Q(N, H_T + _levenbergMarquardt * MayerHessianMatrix::Identity());
			_QP.set_q(N, g_T);
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
