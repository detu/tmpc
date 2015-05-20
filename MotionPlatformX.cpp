#pragma once

#include "MotionPlatformX.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;

static std::vector<MotionLimits> DefaultAxesLimits()
{
	std::vector<MotionLimits> limits;
	limits.push_back(MotionLimits(-1., 1., -1., 1., -1., 1.));
	return limits;
}

MotionPlatformX::MotionPlatformX()
	: MotionPlatform(DefaultAxesLimits())
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
		Eigen::Map<RowMajorMatrix> map_C(C, getOutputDim(), getStateDim());
		map_C.fill(0.);
	}

	// D = dy/du
	if (D)
	{
		Eigen::Map<RowMajorMatrix> map_D(D, getOutputDim(), getInputDim());
		map_D.fill(0.);
		map_D(0, 0) = -1;
	}
}
