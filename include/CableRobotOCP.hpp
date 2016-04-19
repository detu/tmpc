/*
 * CableRobotOCP.h
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#pragma once

#include <OptimalControlProblem.hpp>

namespace mpmc
{
	class CableRobotOCP : public camels::OptimalControlProblem<CableRobotOCP, 13, 8, 13>
	{
	public:
		CableRobotOCP();
		~CableRobotOCP();
		ODEOutput ODE(unsigned t, const StateVector& x, const InputVector& u, const ParamVector& p);

	private:
	};
} /* namespace mpmc */

