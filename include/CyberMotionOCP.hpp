/*
 * CyberMotionOCP.h
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#pragma once

#include <OptimalControlProblem.hpp>

#include "CyberMotion.hpp"

namespace mpmc
{
	class CyberMotionOCP : public camels::OptimalControlProblem<CyberMotionOCP, 2 * CyberMotion::numberOfAxes, CyberMotion::numberOfAxes, 0>
	{
	public:
		ODEOutput ODE(unsigned t, const StateVector& x, const InputVector& u, const ParamVector& p);

	private:
		static const auto N = CyberMotion::numberOfAxes;
		typedef Eigen::Matrix<Scalar, N, N> MatrixNN;

		static StateSensitivityMatrix SensX();
		static InputSensitivityMatrix SensU();
	};
} /* namespace mpmc */

