#pragma once

#include "MotionPlatform.hpp"

#include "casadi_interface/GeneratedFunction.hpp"

namespace mpmc
{
	class CyberMotion : public MotionPlatform
	{
	public:
		CyberMotion();
		void getAxesLimits(Eigen::VectorXd& q_min, Eigen::VectorXd& q_max, Eigen::VectorXd& v_min, Eigen::VectorXd& v_max, Eigen::VectorXd& u_min, Eigen::VectorXd& u_max) const override;
		Eigen::VectorXd getDefaultAxesPosition() const override;
		void Output(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& y, Eigen::MatrixXd& C, Eigen::MatrixXd& D) const override;

		static const unsigned numberOfAxes = 8;

	private:
		mutable casadi_interface::GeneratedFunction _ode;
		mutable casadi_interface::GeneratedFunction _output;
	};
}
