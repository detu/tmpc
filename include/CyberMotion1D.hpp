#pragma once

#include "CyberMotion.hpp"

#include <Eigen/Dense>

namespace mpmc
{
	// Implementation of the CMS kinematics model where only axis 0 (linear track) can move.
	class CyberMotion1D : public MotionPlatform
	{
	public:
		typedef Eigen::Matrix<double, 8, 1> FullPositionVector;
		typedef Eigen::Matrix<double, 16, 1> FullStateVector;
		typedef Eigen::Matrix<double, 8, 1> FullInputVector;
		typedef Eigen::Matrix<double, 6, 1> OutputVector;
		typedef Eigen::Matrix<double, 6, 16, Eigen::ColMajor> FullJacobianMatrixC;
		typedef Eigen::Matrix<double, 6, 8, Eigen::ColMajor> FullJacobianMatrixD;
		typedef Eigen::Matrix<double, 6, 2, Eigen::ColMajor> JacobianMatrixC;
		typedef Eigen::Matrix<double, 6, 1, Eigen::ColMajor> JacobianMatrixD;

		CyberMotion1D();
		void getAxesLimits(Eigen::VectorXd& q_min, Eigen::VectorXd& q_max, Eigen::VectorXd& v_min, Eigen::VectorXd& v_max, Eigen::VectorXd& u_min, Eigen::VectorXd& u_max) const override;
		Eigen::VectorXd getDefaultAxesPosition() const override;
		void Output(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& y, Eigen::MatrixXd& C, Eigen::MatrixXd& D) const override;

		FullStateVector getFullState() const;
		void setFullState(mpmc::CyberMotion1D::FullStateVector val);

	private:
		CyberMotion _fullCMS;
		FullStateVector _fullState;
	};
}