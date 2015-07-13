#include <MultiStageQP.hpp>

#include <gtest/gtest.h>

#include <iostream>

TEST(test_1, my_test_1)
{
	camels::MultiStageQP qp(2, 1, 2);
	qp.H(0).setIdentity();
	qp.H(1).setIdentity();
	qp.H(2).setIdentity();

	qp.C(0) << 1, 1, 0.5,
			   0, 1, 1.0;
	qp.C(1) << 1, 1, 0.5,
			   0, 1, 1.0;

	Eigen::MatrixXd M(4, 4);
	Eigen::VectorXd v(2);

	qp.Calculate_M(M);
	//std::cout << M;

	Eigen::MatrixXd M_expected(4, 4);
	M_expected << 1, 1, 0.5,   0,
				  0, 1,   1,   0,
				  1, 2, 1.5, 0.5,
				  0, 1,   1,   1;
	
	EXPECT_TRUE(M_expected == M);

	//qp.PrintQP_C(std::cout);
}