#pragma once

#include "Trajectory.hpp"
#include "OptimalControlProblem.hpp"
#include "matrix.hpp"
#include "../qp/qp.hpp"

#include <stdexcept>
//#include <type_traits>

namespace tmpc
{
	template <typename OCP, typename Integrator_, class QPSolver_>
	class RealtimeIteration
	{
		static auto const NX = OCP::NX;
		static auto const NU = OCP::NU;

	public:
		typedef QPSolver_ QPSolver;
		typedef Integrator_ Integrator;
		typedef Trajectory<NX, NU> WorkingPoint;

		RealtimeIteration(OCP const& ocp, Integrator const& integrator,
				QPSolver& solver, WorkingPoint const& working_point)
		:	_ocp(ocp)
		,	_QP(working_point.nT())
		,	_workingPoint(working_point)
		,	_solution(working_point.nT())
		,	solver_(solver)
		,	_integrator(integrator)
		,	_prepared(true)
		{
			UpdateQP();
		}

		RealtimeIteration(RealtimeIteration const&) = delete;

		// Feed current state x0, get back control input u.
		template <typename StateVector>
		decltype(auto) Feedback(const StateVector& x0)
		{
			if (!_prepared)
				throw std::logic_error("ModelPredictiveController::Feedback(): controller is not prepared.");

			/** embed current initial value */
			auto const w0 = eval(x0 - _workingPoint.get_x(0));
			_QP.set_x_min(0, w0);
			_QP.set_x_max(0, w0);

			/** solve QP */
			solver_.Solve(_QP, _solution);

			// Add QP step to the working point.
			_workingPoint += _solution;
			_prepared = false;

			// Return the calculated control input.
			return _workingPoint.get_u(0);
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
			typedef typename OCP::StateVector StateVector;
			typedef typename OCP::InputVector InputVector;
			typedef typename OCP::StateInputVector StateInputVector;
			typedef typename OCP::LagrangeHessianMatrix LagrangeHessianMatrix;
			typedef typename OCP::MayerHessianMatrix MayerHessianMatrix;

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
				set_H(_QP, i, H_i + _levenbergMarquardt * identity<LagrangeHessianMatrix>());
				set_g(_QP, i, g_i);

				// Bound constraints.
				StateInputVector z_min, z_max;
				z_min << _ocp.getStateMin(), _ocp.getInputMin();
				z_max << _ocp.getStateMax(), _ocp.getInputMax();

				// C = [ssA, ssB];
				// x_{k+1} = C * z_k + c_k
				StateVector x_plus;

				Eigen::Matrix<double, NX, NX> A;
				Eigen::Matrix<double, NX, NU> B;
				_integrator.Integrate(_ocp.getODE(), i * _integrator.timeStep(), _workingPoint.get_x(i), _workingPoint.get_u(i), x_plus, A, B);
				_QP.set_A(i, A);
				_QP.set_B(i, B);

				// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
				// c = f(z_k) - x_{k+1}
				_QP.set_b(i, x_plus - _workingPoint.get_x(i + 1));

				typename OCP::ConstraintJacobianMatrix D;
				typename OCP::ConstraintVector d_min, d_max;
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
			typename OCP::TerminalConstraintJacobianMatrix D;
			typename OCP::TerminalConstraintVector d_min, d_max;
			_ocp.TerminalConstraints(xN, D, d_min, d_max);
			_QP.set_C_end(D);
			_QP.set_d_end_min(d_min);
			_QP.set_d_end_max(d_max);

			_QP.set_x_min(N, _ocp.getTerminalStateMin() - xN);
			_QP.set_x_max(N, _ocp.getTerminalStateMax() - xN);

			// Hessian and gradient of Mayer term.
			typename OCP::MayerHessianMatrix H_T;
			typename OCP::StateVector g_T;
			_ocp.MayerTerm(xN, g_T, H_T);

			// Adding Levenberg-Marquardt term to make H positive-definite.
			_QP.set_Q(N, H_T + _levenbergMarquardt * identity<MayerHessianMatrix>());
			_QP.set_q(N, g_T);
		}

		// Private data members.
		OCP const& _ocp;
		Integrator const& _integrator;
		
		typename QPSolver::Problem _QP;
		typename QPSolver::Solution _solution;
		QPSolver& solver_;

		double _levenbergMarquardt = 0.;

		// Working point (linearization point).
		WorkingPoint _workingPoint;

		// Preparation() sets this flag to true, Feedback() resets it to false.
		bool _prepared;
	};
}
