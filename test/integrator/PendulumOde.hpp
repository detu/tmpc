#pragma once

#include <tmpc/casadi/GeneratedFunction.hpp>

#include <pendulum_ode_generated.h>

#include <blaze/Math.h>

#include <cstddef>


namespace tmpc :: testing
{
    class PendulumOdeBase
	{
	public:
		static unsigned const NX = 2;
		static unsigned const NU = 1;
		static unsigned const NQ = 2;
		static unsigned const NR = 2;

        typedef blaze::StaticVector<double, NX, blaze::columnVector> StateVector;
		typedef blaze::StaticVector<double, NU, blaze::columnVector> InputVector;
		typedef blaze::StaticVector<double, NQ, blaze::columnVector> QuadVector;
		typedef blaze::StaticVector<double, NR, blaze::columnVector> ResVector;
		typedef blaze::StaticMatrix<double, NX, NX, blaze::columnMajor> StateStateMatrix;
		typedef blaze::StaticMatrix<double, NX, NU, blaze::columnMajor> StateInputMatrix;
		typedef blaze::StaticMatrix<double, NU, NU, blaze::columnMajor> InputInputMatrix;
		typedef blaze::StaticMatrix<double, NQ, NX, blaze::columnMajor> QuadStateMatrix;
		typedef blaze::StaticMatrix<double, NQ, NU, blaze::columnMajor> QuadInputMatrix;
		typedef blaze::StaticMatrix<double, NR, NX, blaze::columnMajor> ResStateMatrix;
		typedef blaze::StaticMatrix<double, NR, NU, blaze::columnMajor> ResInputMatrix;

	protected:
		casadi::GeneratedFunction const _ode {pendulum_ode_functions()};
	};


	class PendulumOde 
	: 	public PendulumOdeBase
	{
	public:
		/**
		 * \brief Evaluates ODE.
		 */
		template <typename VT1, bool TF1, typename MT1, bool SO1, typename MT2, bool SO2>
		void operator()(double t, StateVector const& x0, InputVector const& u0,	
			blaze::Vector<VT1, TF1>& xdot, 
			blaze::Matrix<MT1, SO1>& A, 
			blaze::Matrix<MT2, SO2>& B) const
		{
			std::nullptr_t null;
			_ode(std::tie(t, x0, u0), std::tie(xdot, A, B, null, null, null, null, null, null));
		}


		/**
		 * \brief Evaluates ODE and quadrature.
		 */
		void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
			QuadVector& q, QuadStateMatrix& qA, QuadInputMatrix& qB) const
		{
			std::nullptr_t null;
			_ode(std::tie(t, x0, u0), std::tie(xdot, A, B, q, qA, qB, null, null, null));
		}


		/**
		 * \brief Evaluates ODE, quadrature and residuals.
		 */
		void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
			QuadVector& q, QuadStateMatrix& qA, QuadInputMatrix& qB, ResVector& r, ResStateMatrix& rA, ResInputMatrix& rB) const
		{
			std::nullptr_t null;
			_ode(std::tie(t, x0, u0), std::tie(xdot, A, B, q, qA, qB, r, rA, rB));
		}


		void operator()(double t, StateVector const& x0, InputVector const& u0,
				StateVector const& x0_seed, InputVector const& u_seed, StateVector& xdot, StateVector& xdot_sens) const
		{
			std::nullptr_t null;
			static casadi::GeneratedFunction const _ode(pendulum_ode_sens_functions());
			_ode(std::tie(t, x0, u0, x0_seed, u_seed), std::tie(xdot, xdot_sens));
		}

		/**
		 * \brief Evaluates ODE without sensitivities.
		 */
		StateVector operator()(double t, StateVector const& x0, InputVector const& u0) const
		{
			std::nullptr_t null;
			StateVector xdot;
			_ode(std::tie(t, x0, u0), std::tie(xdot, null, null, null, null, null, null, null, null));

			return xdot;
		}
	};


	class PendulumOdeR 
	: 	public PendulumOdeBase
	{
	public:
		/**
		 * \brief Evaluates ODE and residuals.
		 */
		void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
			ResVector& r, ResStateMatrix& rA, ResInputMatrix& rB) const
		{
			std::nullptr_t null;
			_ode(std::tie(t, x0, u0), std::tie(xdot, A, B, null, null, null, r, rA, rB));
		}
	};
}