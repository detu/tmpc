/*
 * RK4.hpp
 *
 *  Created on: May 8, 2016
 *      Author: kotlyar
 */

#ifndef INCLUDE_INTEGRATOR_RK4_HPP_
#define INCLUDE_INTEGRATOR_RK4_HPP_

#include <Eigen/Dense>

namespace camels
{
	template<unsigned N, typename Matrix>
	decltype(auto) topRows(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template topRows<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) topRows(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template topRows<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) bottomRows(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template bottomRows<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) bottomRows(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template bottomRows<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) leftCols(Eigen::MatrixBase<Matrix>& m)
	{
		return m.template leftCols<N>();
	}

	template<unsigned N, typename Matrix>
	decltype(auto) leftCols(Eigen::MatrixBase<Matrix> const& m)
	{
		return m.template leftCols<N>();
	}

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

			auto const x0 = topRows   <NX>(z0);
			auto const u  = bottomRows<NU>(z0);

			// Calculating next state
			StateInputVector z1;
												_ode.ODE(t0,          z0, k1, J1);
			z1 << x0 + k1 * h / 2., u;			_ode.ODE(t0 + h / 2., z1, k2, J2);
			z1 << x0 + k2 * h / 2., u;			_ode.ODE(t0 + h / 2., z1, k3, J3);
			z1 << x0 + k3 * h     , u;			_ode.ODE(t0 + h,      z1, k4, J4);

			x_next = x0 + (k1 + 2. * k2 + 2. * k3 + k4) * h / 6.;

			// Calculating sensitivities
			ODEJacobianMatrix const& J1_bar = J1;
			ODEJacobianMatrix const  J2_bar = J2 + h / 2. * leftCols<NX>(J2) * J1_bar;
			ODEJacobianMatrix const  J3_bar = J3 + h / 2. * leftCols<NX>(J3) * J2_bar;
			ODEJacobianMatrix const  J4_bar = J4 + h      * leftCols<NX>(J4) * J3_bar;

			J = ODEJacobianMatrix::Identity() + h / 6. * (J1_bar + 2. * J2_bar + 2. * J3_bar + J4_bar);
		}

		double timeStep() const noexcept { return _timeStep; }

	private:
		ODEModel const& _ode;
		double const _timeStep;
	};
}

#endif /* INCLUDE_INTEGRATOR_RK4_HPP_ */
