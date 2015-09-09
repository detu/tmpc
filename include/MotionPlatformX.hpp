#pragma once

#include "MotionPlatform.hpp"

namespace mpmc
{
	class MotionPlatformX : public MotionPlatform
	{
	public:
		MotionPlatformX();
		void getAxesLimits(Eigen::VectorXd& q_min, Eigen::VectorXd& q_max, Eigen::VectorXd& v_min, Eigen::VectorXd& v_max, Eigen::VectorXd& u_min, Eigen::VectorXd& u_max) const override;
		Eigen::VectorXd getDefaultAxesPosition() const override;
		void Output(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& y, Eigen::MatrixXd& C, Eigen::MatrixXd& D) const override;
	};
}