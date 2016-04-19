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
	CableRobotOCP::CableRobotOCP()
	{
		cablerobot_ode_incref();
		cablerobot_output_incref();
	}

	CableRobotOCP::~CableRobotOCP()
	{
		cablerobot_ode_decref();
		cablerobot_output_decref();
	}

	/*
	CableRobotOCP::ODEOutput CableRobotOCP::ODE(unsigned t, const StateVector& x, const InputVector& u, const ParamVector& p)
	{
		StateVector xdot;
		xdot << x.bottomRows(N), u;

		static const StateSensitivityMatrix sens_x = SensX();
		static const InputSensitivityMatrix sens_u = SensU();

		return ODEOutput(xdot, sens_x, sens_u);
	}
	*/
} /* namespace mpmc */
