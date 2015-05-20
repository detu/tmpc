#pragma once

#include "MotionPlatform.hpp"

class MotionPlatformX : public MotionPlatform
{
public:
	MotionPlatformX();
	void Output(const ConstVectorMap& x, const ConstVectorMap& u, VectorMap& y, MatrixMap& C, MatrixMap& D) const override;
};