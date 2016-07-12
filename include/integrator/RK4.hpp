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
	class RK4
	{
	public:
		RK4(double time_step) : _timeStep(time_step) {}

		template <typename ODE, typename StateVector0_, typename InputVector_, typename StateVector1_, typename AMatrix, typename BMatrix>
		void Integrate(ODE const& ode, double t0, StateVector0_ const& x0, InputVector_ const& u, StateVector1_& x_next, AMatrix& A, BMatrix& B) const
		{
			auto constexpr NX = rows<StateVector0_>();
			auto constexpr NU = rows<InputVector_ >();

			typedef Eigen::Matrix<double, NX,  1> StateVector;
			typedef Eigen::Matrix<double, NX, NX> StateStateMatrix;
			typedef Eigen::Matrix<double, NX, NU> StateInputMatrix;

			StateVector k1, k2, k3, k4;
			StateStateMatrix A1, A2, A3, A4;
			StateInputMatrix B1, B2, B3, B4;
			auto const h = _timeStep;

			// Calculating next state
			ode(t0,          x0              , u, k1, A1, B1);
			ode(t0 + h / 2., x0 + k1 * h / 2., u, k2, A2, B2);
			ode(t0 + h / 2., x0 + k2 * h / 2., u, k3, A3, B3);
			ode(t0 + h,      x0 + k3 * h     , u, k4, A4, B4);

			x_next = x0 + (k1 + 2. * k2 + 2. * k3 + k4) * h / 6.;

			// Calculating sensitivities
			auto const& A1_bar =      A1;							auto const& B1_bar =      B1;
			auto const  A2_bar = eval(A2 + h / 2. * A2 * A1_bar);	auto const  B2_bar = eval(B2 + h / 2. * A2 * B1_bar);
			auto const  A3_bar = eval(A3 + h / 2. * A3 * A2_bar);	auto const  B3_bar = eval(B3 + h / 2. * A3 * B2_bar);
			auto const  A4_bar =      A4 + h      * A4 * A3_bar ;	auto const  B4_bar =      B4 + h      * A4 * B3_bar ;

			A = identity<StateStateMatrix>() + h / 6. * (A1_bar + 2. * A2_bar + 2. * A3_bar + A4_bar);
			B = 					           h / 6. * (B1_bar + 2. * B2_bar + 2. * B3_bar + B4_bar);
		}

		//
		// Sensitivity-free version of the Integrate() function.
		//
		template <typename ODE, typename StateVector, typename InputVector>
		decltype(auto) Integrate(ODE const& ode, double t0, StateVector const& x0, InputVector const& u) const
		{
			auto const h = _timeStep;

			// Calculating next state
			auto const k1 = eval(ode(t0,          x0              , u));
			auto const k2 = eval(ode(t0 + h / 2., x0 + k1 * h / 2., u));
			auto const k3 = eval(ode(t0 + h / 2., x0 + k2 * h / 2., u));
			auto const k4 = eval(ode(t0 + h,      x0 + k3 * h     , u));

			return x0 + (k1 + 2. * k2 + 2. * k3 + k4) * h / 6.;
		}

		// TODO: should it be a parameter of Integrate() instead?
		double timeStep() const noexcept { return _timeStep; }

	private:
		double const _timeStep;
	};
}

#endif /* INCLUDE_INTEGRATOR_RK4_HPP_ */
