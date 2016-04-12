/*
 * CyberMotionOCP.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#include <CyberMotionOCP.hpp>

namespace mpmc
{
	CyberMotionOCP::ODEOutput CyberMotionOCP::ODE(const StateVector& x, const InputVector& u)
	{
		StateVector xdot;
		xdot << x.bottomRows(N), u;

		static const StateSensitivityMatrix sens_x = SensX();
		static const InputSensitivityMatrix sens_u = SensU();

		return ODEOutput(xdot, sens_x, sens_u);
	}

	CyberMotionOCP::StateSensitivityMatrix CyberMotionOCP::SensX()
	{
		StateSensitivityMatrix sens_x;
		sens_x << MatrixNN::Zero(), MatrixNN::Identity(), MatrixNN::Zero(), MatrixNN::Zero();
		return sens_x;
	}

	CyberMotionOCP::InputSensitivityMatrix CyberMotionOCP::SensU()
	{
		InputSensitivityMatrix sens_u;
		sens_u << MatrixNN::Zero(), MatrixNN::Identity();
		return sens_u;
	}
} /* namespace mpmc */
