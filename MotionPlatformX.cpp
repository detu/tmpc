#pragma once

#include "MotionPlatformX.hpp"

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

void MotionPlatformX::Output(const ConstVectorMap& x, const ConstVectorMap& u, VectorMap& y, MatrixMap& C, MatrixMap& D) const
{
	if (!(x.size() == getStateDim() && u.size() == getInputDim() && y.size() == getOutputDim()
		&& C.rows() == getOutputDim() && C.cols() == getStateDim() && D.rows() == getOutputDim() && D.cols() == getInputDim()))
		throw std::invalid_argument("MotionPlatformX::Output() input argument dimension mismatch");

	// Specific force
	y.block<3, 1>(0, 0) = getGravity();
	y(0) -= u(0);

	// Rotational velocity
	y.block<3, 1>(3, 0).fill(0.);

	// C = dy/dx
	C.fill(0.);

	// D = dy/du
	D.fill(0.);
	D(0, 0) = -1;
}
