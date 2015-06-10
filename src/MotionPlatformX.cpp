#pragma once

#include "MotionPlatformX.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ColMajorMatrix;

namespace rtmc
{
	MotionPlatformX::MotionPlatformX()
		: MotionPlatform(1)
	{
	}

	void MotionPlatformX::Output(const double * x, const double * u, double * y, double * C, double * D) const
	{
		// Specific force
		Eigen::Map<Eigen::Vector3d> f(y);
		f = getGravity();
		f(0) -= u[0];

		// Rotational velocity
		Eigen::Map<Eigen::Vector3d> omega(y + 3);
		omega.fill(0.);

		// C = dy/dx
		if (C)
		{
			Eigen::Map<ColMajorMatrix> map_C(C, getOutputDim(), getStateDim());
			map_C.fill(0.);
		}

		// D = dy/du
		if (D)
		{
			Eigen::Map<ColMajorMatrix> map_D(D, getOutputDim(), getInputDim());
			map_D.fill(0.);
			map_D(0, 0) = -1;
		}
	}

	void MotionPlatformX::getDefaultAxesPosition(double * q) const
	{
		q[0] = 0;
	}

	void MotionPlatformX::getAxesLimits(double * q_min, double * q_max, double * v_min, double * v_max, double * u_min, double * u_max) const
	{
		q_min[0] = -1;	q_max[0] = 1;
		v_min[0] = -1;	v_max[0] = 1;
		u_min[0] = -1;  u_max[0] = 1;
	}
}