/*
 * CableRobotOCP.h
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#pragma once

#include <OptimalControlProblem.hpp>

#include "CasADiGeneratedFunction.hpp"

namespace mpmc
{
	class CableRobotOCP : public camels::OptimalControlProblem<CableRobotOCP, 13, 8>
	{
	public:
		CableRobotOCP(unsigned Nt);

		void ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac);

	private:
		mutable CasADiGeneratedFunction _ode;
		mutable CasADiGeneratedFunction _output;
	};
} /* namespace mpmc */

