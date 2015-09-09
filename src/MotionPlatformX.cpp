#include "MotionPlatformX.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ColMajorMatrix;

namespace mpmc
{
	MotionPlatformX::MotionPlatformX()
		: MotionPlatform(1)
	{
	}

	void MotionPlatformX::Output(const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::VectorXd& y, Eigen::MatrixXd& C, Eigen::MatrixXd& D) const
	{
		// Specific force
		auto f = getGravity();
		f(0) -= u[0];

		// Rotational velocity
		static const Eigen::Vector3d omega(0., 0., 0.);

		y << f, omega;

		// C = dy/dx
		C << Eigen::MatrixXd::Zero(getOutputDim(), getStateDim());

		// D = dy/du
		D << Eigen::MatrixXd::Zero(getOutputDim(), getInputDim());
		D(0, 0) = -1;
	}

	Eigen::VectorXd MotionPlatformX::getDefaultAxesPosition() const
	{
		Eigen::VectorXd result(1);
		result << 0.;
		return result;
	}

	void MotionPlatformX::getAxesLimits(Eigen::VectorXd& q_min, Eigen::VectorXd& q_max, Eigen::VectorXd& v_min, Eigen::VectorXd& v_max, Eigen::VectorXd& u_min, Eigen::VectorXd& u_max) const
	{
		q_min << -1;	q_max << 1;
		v_min << -1;	v_max << 1;
		u_min << -1;	u_max << 1;
	}
}
