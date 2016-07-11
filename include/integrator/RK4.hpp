/*
 * RK4.hpp
 *
 *  Created on: May 8, 2016
 *      Author: kotlyar
 */

#ifndef INCLUDE_INTEGRATOR_RK4_HPP_
#define INCLUDE_INTEGRATOR_RK4_HPP_

#include <Eigen/Dense>

#include "../core/matrix.hpp"

namespace tmpc
{
	template<typename ODEModel_>
	class RK4
	{
	public:
		typedef ODEModel_ ODEModel;
		typedef typename ODEModel::StateVector StateVector;
		typedef typename ODEModel::StateInputVector StateInputVector;
		typedef typename ODEModel::ODEJacobianMatrix ODEJacobianMatrix;

		static unsigned const NX = ODEModel::NX;
		static unsigned const NU = ODEModel::NU;

		RK4(ODEModel const& ode, double time_step) : _ode(ode), _timeStep(time_step) {}

		void Integrate(double t0, StateInputVector const& z0, StateVector& x_next, ODEJacobianMatrix& J) const
		{
			StateVector k1, k2, k3, k4;
			ODEJacobianMatrix J1, J2, J3, J4;
			auto const h = _timeStep;

			auto const x0 = top_rows   <NX>(z0);
			auto const u  = bottom_rows<NU>(z0);

			// Calculating next state
			StateInputVector z1;
												_ode.ODE(t0,          z0, k1, J1);
			z1 << x0 + k1 * h / 2., u;			_ode.ODE(t0 + h / 2., z1, k2, J2);
			z1 << x0 + k2 * h / 2., u;			_ode.ODE(t0 + h / 2., z1, k3, J3);
			z1 << x0 + k3 * h     , u;			_ode.ODE(t0 + h,      z1, k4, J4);

			x_next = x0 + (k1 + 2. * k2 + 2. * k3 + k4) * h / 6.;

			// Calculating sensitivities
			ODEJacobianMatrix const& J1_bar = J1;
			ODEJacobianMatrix const  J2_bar = J2 + h / 2. * left_cols<NX>(J2) * J1_bar;
			ODEJacobianMatrix const  J3_bar = J3 + h / 2. * left_cols<NX>(J3) * J2_bar;
			ODEJacobianMatrix const  J4_bar = J4 + h      * left_cols<NX>(J4) * J3_bar;

			J = ODEJacobianMatrix::Identity() + h / 6. * (J1_bar + 2. * J2_bar + 2. * J3_bar + J4_bar);
		}

		template <typename StateVector0_, typename InputVector_, typename StateVector1_, typename AMatrix, typename BMatrix>
		void Integrate(double t0, StateVector0_ const& x0, InputVector_ const& u, StateVector1_& x_next, AMatrix& A, BMatrix& B) const
		{
			typedef Eigen::Matrix<double, NX,  1> StateVector;
			typedef Eigen::Matrix<double, NX, NX> StateStateMatrix;
			typedef Eigen::Matrix<double, NX, NU> StateInputMatrix;

			StateVector k1, k2, k3, k4;
			StateStateMatrix A1, A2, A3, A4;
			StateInputMatrix B1, B2, B3, B4;
			auto const h = _timeStep;

			// Calculating next state
			_ode.ODE(t0,          x0              , u, k1, A1, B1);
			_ode.ODE(t0 + h / 2., x0 + k1 * h / 2., u, k2, A2, B2);
			_ode.ODE(t0 + h / 2., x0 + k2 * h / 2., u, k3, A3, B3);
			_ode.ODE(t0 + h,      x0 + k3 * h     , u, k4, A4, B4);

			x_next = x0 + (k1 + 2. * k2 + 2. * k3 + k4) * h / 6.;

			// Calculating sensitivities
			StateStateMatrix const& A1_bar = A1;						StateInputMatrix const& B1_bar = B1;
			StateStateMatrix const  A2_bar = A2 + h / 2. * A2 * A1_bar;	StateInputMatrix const  B2_bar = B2 + h / 2. * A2 * B1_bar;
			StateStateMatrix const  A3_bar = A3 + h / 2. * A3 * A2_bar;	StateInputMatrix const  B3_bar = B3 + h / 2. * A3 * B2_bar;
			StateStateMatrix const  A4_bar = A4 + h      * A4 * A3_bar;	StateInputMatrix const  B4_bar = B4 + h      * A4 * B3_bar;

			A = identity<StateStateMatrix>() + h / 6. * (A1_bar + 2. * A2_bar + 2. * A3_bar + A4_bar);
			B = 					           h / 6. * (B1_bar + 2. * B2_bar + 2. * B3_bar + B4_bar);
		}

		double timeStep() const noexcept { return _timeStep; }

	private:
		ODEModel const& _ode;
		double const _timeStep;
	};
}

#endif /* INCLUDE_INTEGRATOR_RK4_HPP_ */
