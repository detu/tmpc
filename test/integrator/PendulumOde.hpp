#pragma once

#include <tmpc/casadi_interface/GeneratedFunction.hpp>

#include <pendulum_ode_generated.h>

#include <blaze/Math.h>


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
		casadi_interface::GeneratedFunction const _ode {pendulum_ode_functions()};
	};


	class PendulumOde 
	: 	public PendulumOdeBase
	{
	public:
		/**
		 * \brief Evaluates ODE.
		 */
		template <typename VT1, typename MT1, typename MT2>
		void operator()(double t, StateVector const& x0, InputVector const& u0,	
			blaze::DenseVector<VT1, blaze::columnVector>& xdot, 
			blaze::DenseMatrix<MT1, blaze::columnMajor>& A, 
			blaze::DenseMatrix<MT2, blaze::columnMajor>& B) const
		{
			_ode({&t, x0.data(), u0.data()}, {data(xdot), data(A), data(B), nullptr, nullptr, nullptr, nullptr, nullptr, nullptr});
		}


		/**
		 * \brief Evaluates ODE.
		 */
		template <typename VT1, bool TF1, typename MT1, bool SO1, typename MT2, bool SO2>
		void operator()(double t, StateVector const& x0, InputVector const& u0,	
			blaze::Vector<VT1, TF1>& xdot, 
			blaze::Matrix<MT1, SO1>& A, 
			blaze::Matrix<MT2, SO2>& B) const
		{
			_ode({&t, x0.data(), u0.data()}, {xdot_.data(), A_.data(), B_.data(), nullptr, nullptr, nullptr, nullptr, nullptr, nullptr});
			
			~xdot = xdot_;
			~A = A_;
			~B = B_;
		}


		/**
		 * \brief Evaluates ODE and quadrature.
		 */
		void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
			QuadVector& q, QuadStateMatrix& qA, QuadInputMatrix& qB) const
		{
			_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data(), q.data(), qA.data(), qB.data(), nullptr, nullptr, nullptr});
		}

		/**
		 * \brief Evaluates ODE, quadrature and residuals.
		 */
		void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
			QuadVector& q, QuadStateMatrix& qA, QuadInputMatrix& qB, ResVector& r, ResStateMatrix& rA, ResInputMatrix& rB) const
		{
			_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data(), q.data(), qA.data(), qB.data(), r.data(), rA.data(), rB.data()});
		}

		void operator()(double t, StateVector const& x0, InputVector const& u0,
				StateVector const& x0_seed, InputVector const& u_seed, StateVector& xdot, StateVector& xdot_sens) const
		{
			static casadi_interface::GeneratedFunction const _ode(pendulum_ode_sens_functions());
			_ode({&t, x0.data(), u0.data(), x0_seed.data(), u_seed.data()}, {xdot.data(), xdot_sens.data()});
		}

		/**
		 * \brief Evaluates ODE without sensitivities.
		 */
		StateVector operator()(double t, StateVector const& x0, InputVector const& u0) const
		{
			StateVector xdot;
			_ode({&t, x0.data(), u0.data()}, {xdot.data(), nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr});

			return xdot;
		}


	private:
		mutable StateVector xdot_;
		mutable StateStateMatrix A_;
		mutable StateInputMatrix B_;
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
			_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data(), nullptr, nullptr, nullptr, r.data(), rA.data(), rB.data()});
		}
	};
}