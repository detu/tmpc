/*
 * ODEModel.hpp
 *
 *  Created on: May 4, 2016
 *      Author: kotlyar
 */

#ifndef INCLUDE_CORE_ODEMODEL_HPP_
#define INCLUDE_CORE_ODEMODEL_HPP_

namespace camels
{
	template<typename Derived_>
	class ODEModel
	{
	public:
		typedef Derived_ Derived;

		template<typename StateInputVector, typename StateVector, typename ODEJacobianMatrix>
		void ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac) const
		{
			derived().ODE(t, z, xdot, jac);
		}

	private:
		Derived& derived() { return static_cast<Derived&>(*this); }
		Derived const& derived() const { return static_cast<Derived const&>(*this); }
	};
}

#endif /* INCLUDE_CORE_ODEMODEL_HPP_ */
