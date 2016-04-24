/*
 * CableRobotOCP.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#include <CableRobotOCP.hpp>

#include "cablerobot_generated.h"

namespace mpmc
{
	CableRobotOCP::CableRobotOCP(unsigned Nt)
	:	camels::OptimalControlProblem<CableRobotOCP, 13, 8>(Nt)
	,	_ode(CASADI_GENERATED_FUNCTION_INTERFACE(cablerobot_ode))
	,	_output(CASADI_GENERATED_FUNCTION_INTERFACE(cablerobot_output))
	{

	}

	void CableRobotOCP::ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac)
	{
		_ode({z.data(), nullptr}, {xdot.data(), jac.data()});
	}
} /* namespace mpmc */
