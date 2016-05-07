/*
 * CableRobotOCP.h
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#pragma once

#include <core/OptimalControlProblem.hpp>

#include "casadi_interface/GeneratedFunction.hpp"

namespace mpmc
{
	class CableRobotOCP : public camels::OptimalControlProblem<CableRobotOCP, 13, 8>
	{
	public:
		CableRobotOCP(unsigned Nt);

		void ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac);

	private:
		casadi_interface::GeneratedFunction mutable _ode;
		casadi_interface::GeneratedFunction mutable _output;
	};
} /* namespace mpmc */

