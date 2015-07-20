#include <MultiStageQP.hpp>

#include <gtest/gtest.h>

#include <iostream>

TEST(test_1, my_test_1)
{
	camels::MultiStageQP qp(2, 1, 2);

	// Stage 0
	Eigen::MatrixXd H0(qp.nZ(), qp.nZ());
	H0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H0 = H0.transpose() * H0;	// Make positive definite.

	const Eigen::MatrixXd Q0 = H0.topLeftCorner(qp.nX(), qp.nX());
	const Eigen::MatrixXd R0 = H0.bottomRightCorner(qp.nU(), qp.nU());
	const Eigen::MatrixXd S0 = H0.topRightCorner(qp.nX(), qp.nU());
	const Eigen::MatrixXd S0T = H0.bottomLeftCorner(qp.nU(), qp.nX());

	Eigen::MatrixXd A0(qp.nX(), qp.nX());
	A0 << 1, 1, 0, 1;

	Eigen::MatrixXd B0(qp.nX(), qp.nU());
	B0 << 0.5, 1.0;

	Eigen::VectorXd a0(qp.nX());
	a0 << 1, 2;

	// Stage 1
	Eigen::MatrixXd H1(qp.nZ(), qp.nZ());
	H1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H1 = H1.transpose() * H1;	// Make positive definite.

	const Eigen::MatrixXd Q1 = H1.topLeftCorner(qp.nX(), qp.nX());
	const Eigen::MatrixXd R1 = H1.bottomRightCorner(qp.nU(), qp.nU());
	const Eigen::MatrixXd S1 = H1.topRightCorner(qp.nX(), qp.nU());
	const Eigen::MatrixXd S1T = H1.bottomLeftCorner(qp.nU(), qp.nX());

	Eigen::MatrixXd A1(qp.nX(), qp.nX());
	A1 << 1, 1, 0, 1;

	Eigen::MatrixXd B1(qp.nX(), qp.nU());
	B1 << 0.5, 1.0;

	Eigen::VectorXd a1(qp.nX());
	a1 << 1, 2;

	// Stage 2
	Eigen::MatrixXd H2(qp.nX(), qp.nX());
	H2 << 1, 2, 3, 4;
	H2 = H2.transpose() * H2;	// Make positive definite.

	const Eigen::MatrixXd Q2 = H2.topLeftCorner(qp.nX(), qp.nX());

	// Setup QP
	qp.H(0) = H0;
	qp.H(1) = H1;
	qp.H(2) = H2;

	qp.C(0) << A0, B0;
	qp.c(0) << a0;
	qp.C(1) << A1, B1;
	qp.c(1) << a1;

	// Condense
	Eigen::MatrixXd M(qp.nDep(), qp.nIndep());
	Eigen::VectorXd v(qp.nDep());
	Eigen::SparseMatrix<double> H(qp.nVar(), qp.nVar());
	Eigen::SparseMatrix<double> P(qp.nVar(), qp.nVar());
	Eigen::MatrixXd Hc(qp.nIndep(), qp.nIndep());
	Eigen::VectorXd gc(qp.nIndep());

	qp.Calculate_M(M);
	qp.Calculate_v(v);
	qp.Condense(Hc, gc);
// 	std::cout << M << std::endl;
// 	std::cout << v << std::endl;
// 	std::cout << H << std::endl;
// 	std::cout << P << std::endl;
// 	std::cout << Hc << std::endl;
// 	std::cout << gc << std::endl;

	Eigen::MatrixXd M_expected(qp.nDep(), qp.nIndep());
	M_expected << 
		A0,			B0,			Eigen::MatrixXd::Zero(qp.nX(), qp.nU()),
		A1 * A0,	A1 * B0,	B1;

	Eigen::VectorXd v_expected(qp.nIndep());
	v_expected << 
		a0,
		A1 * a0 + a1;
	
	EXPECT_TRUE(M_expected == M);
	EXPECT_TRUE(v_expected == v);

	Eigen::MatrixXd Hc_expected(qp.nIndep(), qp.nIndep());
	
	Hc_expected <<
		A0.transpose() * Q1 * A0 + A0.transpose() * A1.transpose() * Q2 * A1 * A0 + Q0,				A0.transpose() * Q1 * B0 + A0.transpose() * A1.transpose() * Q2 * A1 * B0 + S0,	A0.transpose() * S1 + A0.transpose() * A1.transpose() * Q2 * B1,
		B0.transpose() * Q1 * A0 + B0.transpose() * A1.transpose() * Q2 * A1 * A0 + S0.transpose(), B0.transpose() * Q1 * B0 + B0.transpose() * A1.transpose() * Q2 * A1 * B0 + R0, B0.transpose() * S1 + B0.transpose() * A1.transpose() * Q2 * B1,
		S1.transpose() * A0 + B1.transpose() * Q2 * A1 * A0,										S1.transpose() * B0 + B1.transpose() * Q2 * A1 * B0,							B1.transpose() * Q2 * B1 + R1;

	EXPECT_TRUE(Hc_expected == Hc);
	//std::cout << Hc_expected << std::endl;
	//qp.PrintQP_C(std::cout);

	Eigen::VectorXd gc_expected(qp.nIndep());
	gc_expected <<
		2 * A0.transpose() * Q1 * a0					+ 2 * A0.transpose() * A1.transpose() * Q2 * a1			+ 2 * A0.transpose() * A1.transpose() * Q2 * A1 * a0,
		2 * B0.transpose() * A1.transpose() * Q2 * a1	+ 2 * B0.transpose() * A1.transpose() * Q2 * A1 * a0	+ 2 * B0.transpose() * Q1 * a0,
		2 * B1.transpose() * Q2 * A1 * a0				+ 2 * B1.transpose() * Q2 * a1							+ 2 * S1.transpose() * a0;

	EXPECT_TRUE(gc_expected == gc);
 }