#pragma once

#include "Trajectory.hpp"
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
		static auto const NC = OCP::NC;
		static auto const NCT = OCP::NCT;

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
			auto const N = _ocp.getNumberOfIntervals();

			for (unsigned i = 0; i < N; ++i)
			{
				// Hessian and gradient of Lagrange term.
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
				_ocp.LagrangeTerm(i, _workingPoint.get_x(i), _workingPoint.get_u(i), Q, R, S, q, r);

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
				integrate(_integrator, _ocp.getODE(), i * _integrator.timeStep(), _workingPoint.get_x(i), _workingPoint.get_u(i), x_plus, A, B);
				_QP.set_A(i, A);
				_QP.set_B(i, B);

				// \Delta x_{k+1} = C \Delta z_k + f(z_k) - x_{k+1}
				// c = f(z_k) - x_{k+1}
				_QP.set_b(i, x_plus - _workingPoint.get_x(i + 1));

				Eigen::Matrix<double, NC, NX + NU> D;
				Eigen::Matrix<double, NC, 1> d_min, d_max;
				
				{
					Eigen::Matrix<double, NX + NU, 1> z_i;
					z_i << _workingPoint.get_x(i), _workingPoint.get_u(i);
					_ocp.PathConstraints(i, z_i, D, d_min, d_max);
				}

				set_CD(_QP, i, D);
				_QP.set_d_min(i, d_min);
				_QP.set_d_max(i, d_max);
				
				// Bound constraints.
				_QP.set_x_min(i, _ocp.getStateMin() - _workingPoint.get_x(i));	_QP.set_u_min(i, _ocp.getInputMin() - _workingPoint.get_u(i));
				_QP.set_x_max(i, _ocp.getStateMax() - _workingPoint.get_x(i));	_QP.set_u_max(i, _ocp.getInputMax() - _workingPoint.get_u(i));
			}

			auto const xN = _workingPoint.get_x(N);

			// End state constraints.
			Eigen::Matrix<double, NCT, NX> D;
			Eigen::Matrix<double, NCT, 1> d_min, d_max;
			_ocp.TerminalConstraints(xN, D, d_min, d_max);
			_QP.set_C_end(D);
			_QP.set_d_end_min(d_min);
			_QP.set_d_end_max(d_max);

			_QP.set_x_min(N, _ocp.getTerminalStateMin() - xN);
			_QP.set_x_max(N, _ocp.getTerminalStateMax() - xN);

			// Hessian and gradient of Mayer term.
			Eigen::Matrix<double, NX, NX> H_T;
			Eigen::Matrix<double, NX, 1> g_T;
			_ocp.MayerTerm(xN, g_T, H_T);

			// Adding Levenberg-Marquardt term to make H positive-definite.
			_QP.set_Q(N, H_T + _levenbergMarquardt * identity<Eigen::Matrix<double, NX, NX>>());
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
