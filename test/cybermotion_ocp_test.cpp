#include <CyberMotionOCP.hpp>

#include <gtest/gtest.h>
//#define EXPECT_TRUE(X) assert(X)

#include <iostream>

TEST(test_1, cybermotion_ocp_test)
{
	using mpmc::CyberMotionOCP;
	CyberMotionOCP ocp(1);

	CyberMotionOCP::InputVector u;
	u << 1, 2, 3, 4, 5, 6, 7, 8;

	CyberMotionOCP::StateVector x;
	x << 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24;

	CyberMotionOCP::StateInputVector z;
	z << x, u;

	CyberMotionOCP::StateVector xdot;
	CyberMotionOCP::ODEJacobianMatrix jac;

	ocp.ODE(0, z, xdot, jac);

	const auto N = mpmc::CyberMotion::numberOfAxes;
	typedef Eigen::Matrix<CyberMotionOCP::Scalar, N, N> MatrixNN;

	EXPECT_TRUE(xdot.topRows(N) == x.bottomRows(N));
	EXPECT_TRUE(xdot.bottomRows(N) == u);

	auto const A = jac.leftCols(2 * N);
	auto const B = jac.rightCols(N);

	EXPECT_TRUE(A.topLeftCorner    (N, N) == MatrixNN::Zero());
	EXPECT_TRUE(A.topRightCorner   (N, N) == MatrixNN::Identity());
	EXPECT_TRUE(A.bottomLeftCorner (N, N) == MatrixNN::Zero());
	EXPECT_TRUE(A.bottomRightCorner(N, N) == MatrixNN::Zero());
	EXPECT_TRUE(B.topRows   (N) == MatrixNN::Zero());
	EXPECT_TRUE(B.bottomRows(N) == MatrixNN::Identity());
}
