/*
 * CyberMotionOCP.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#include <CyberMotionOCP.hpp>

#include "cybermotion_generated.h"

namespace mpmc
{
	CyberMotionOCP::CyberMotionOCP()
	:	_ode(CASADI_GENERATED_FUNCTION_INTERFACE(cybermotion_ode))
	,	_output(CASADI_GENERATED_FUNCTION_INTERFACE(cybermotion_output))
	{
		// Initialize state and input limits.

		static const double limits[NU][6] = {
			{ 0.2, 9.3, -1.47, 1.47, -1.0780, 1.0780 },
			{ -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), -1.1802, 1.1802, -1.6762, 1.6762 },
			{ -2.1447, -0.8727, -0.9749, 0.9749, -1.1973, 1.1973 },
			{ -0.7540, 1.5415, -1.1802, 1.1802, -2.1893, 2.1893 },
			{ -3.0159, 3.0159, -1.2999, 1.2999, -0.5644, 0.5644 },
			{ -0.8727, 0.8727, -1.2999, 1.2999, -1.6249, 1.6249 },
			{ -3.0159, 3.0159, -2.0525, 2.0525, -1.3170, 1.3170 },
			// Full-range position for cabin axis is { 0.1938, 1.4839, -0.2450, 0.2450, -0.9800, 0.9800 }.
			// Below, the cabin axis position is limited to both rollers on the curve.
			// This is the only case implemented in current CasADi-generated code.
			{ 0.5700, 1.1603, -0.2450, 0.2450, -0.9800, 0.9800 }
		};

		for (int i = 0; i < NU; ++i)
		{
			_x_min[i] 	   = limits[i][0];	_x_max[i] 	   = limits[i][1];
			_x_min[NU + i] = limits[i][2];	_x_max[NU + i] = limits[i][3];
			_u_min[i] 	   = limits[i][4];	_u_max[i]	   = limits[i][5];
		}

		// The "agile" position and zero velocity.
		_x_default << 4.8078, 0.1218, -1.5319, 0.4760, 0.0006, 0.1396, -0.0005, 0.7991, 0., 0., 0., 0., 0., 0., 0., 0.;
	}

	void CyberMotionOCP::ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac) const
	{
		_ode({z.data(), nullptr}, {xdot.data(), jac.data()});
	}

	void CyberMotionOCP::Output(unsigned t, StateInputVector const& z, OutputVector& y, OutputJacobianMatrix& jac) const
	{
		_output({z.data(), nullptr}, {y.data(), jac.data()});
	}

	CyberMotionOCP::StateVector const& CyberMotionOCP::getStateMin() const
	{
		return _x_min;
	}

	CyberMotionOCP::StateVector const& CyberMotionOCP::getStateMax() const
	{
		return _x_max;
	}

	CyberMotionOCP::InputVector const& CyberMotionOCP::getInputMin() const
	{
		return _u_min;
	}

	CyberMotionOCP::InputVector const& CyberMotionOCP::getInputMax() const
	{
		return _u_max;
	}

	CyberMotionOCP::StateVector const& CyberMotionOCP::getDefaultState() const
	{
		return _x_default;
	}

	//unsigned const CyberMotionOCP::NY;
} /* namespace mpmc */
