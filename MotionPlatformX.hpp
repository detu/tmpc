#pragma once

#include "MotionPlatform.hpp"

class MotionPlatformX : public MotionPlatform
{
public:
	MotionPlatformX();
	void Output(const double * x, const double * u, double * y, double * C = nullptr, double * D = nullptr) const override;
};