#include <CyberMotionOCP.hpp>

#include <gtest/gtest.h>
//#define EXPECT_TRUE(X) assert(X)

#include <iostream>

TEST(test_1, cybermotion_ocp_test)
{
	using mpmc::CyberMotionOCP;
	CyberMotionOCP ocp;

	CyberMotionOCP::InputVector u;
	u << 1, 2, 3, 4, 5, 6, 7, 8;

	CyberMotionOCP::StateVector x;
	x << 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24;

	CyberMotionOCP::ParamVector p;

	const auto ode_out = ocp.ODE(0, x, u, p);

	const auto N = mpmc::CyberMotion::numberOfAxes;
	typedef Eigen::Matrix<CyberMotionOCP::Scalar, N, N> MatrixNN;

	EXPECT_TRUE(ode_out.xDot().topRows(N) == x.bottomRows(N));
	EXPECT_TRUE(ode_out.xDot().bottomRows(N) == u);
	EXPECT_TRUE(ode_out.sensX().topLeftCorner    (N, N) == MatrixNN::Zero());
	EXPECT_TRUE(ode_out.sensX().topRightCorner   (N, N) == MatrixNN::Identity());
	EXPECT_TRUE(ode_out.sensX().bottomLeftCorner (N, N) == MatrixNN::Zero());
	EXPECT_TRUE(ode_out.sensX().bottomRightCorner(N, N) == MatrixNN::Zero());
	EXPECT_TRUE(ode_out.sensU().topRows   (N) == MatrixNN::Zero());
	EXPECT_TRUE(ode_out.sensU().bottomRows(N) == MatrixNN::Identity());
}
