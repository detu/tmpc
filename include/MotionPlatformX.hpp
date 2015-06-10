#pragma once

#include "MotionPlatform.hpp"

namespace rtmc
{
	class MotionPlatformX : public MotionPlatform
	{
	public:
		MotionPlatformX();
		void getAxesLimits(double * q_min, double * q_max, double * v_min, double * v_max, double * u_min, double * u_max) const override;
		void getDefaultAxesPosition(double * q) const override;
		void Output(const double * x, const double * u, double * y, double * C = nullptr, double * D = nullptr) const override;
	};
}