#include <MultiStageQP.hpp>

#include <gtest/gtest.h>

#include <iostream>

TEST(test_1, my_test_1)
{
	camels::MultiStageQP qp(2, 1, 2);
	qp.H(0) << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	qp.H(0) = qp.H(0).transpose() * qp.H(0);
	qp.H(1) << 1, 3, 2, 4, 6, 5, 8, 7, 9;
	qp.H(1) = qp.H(1).transpose() * qp.H(1);
	qp.H(2) << 1, 2, 3, 4;
	qp.H(2) = qp.H(2).transpose() * qp.H(2);

	qp.C(0) << 1, 1, 0.5,
			   0, 1, 1.0;
	qp.c(0) << 1, 2;
	qp.C(1) << 1, 1, 0.5,
			   0, 1, 1.0;
	qp.c(1) << 1, 2;

	Eigen::MatrixXd M(qp.nDep(), qp.nIndep());
	Eigen::VectorXd v(qp.nDep());
	Eigen::SparseMatrix<double> H(qp.nVar(), qp.nVar());
	Eigen::SparseMatrix<double> P(qp.nVar(), qp.nVar());
	Eigen::MatrixXd Hc(qp.nIndep(), qp.nIndep());
	Eigen::VectorXd gc(qp.nIndep());

	qp.Calculate_M(M);
	qp.Calculate_v(v);
	qp.Calculate_FullH(H);
	qp.Calculate_PermutationMatrix(P);
	qp.Condense(Hc, gc);
	std::cout << M << std::endl;
	std::cout << v << std::endl;
	std::cout << H << std::endl;
	std::cout << P << std::endl;
	std::cout << Hc << std::endl;

	Eigen::MatrixXd M_expected(4, 4);
	M_expected << 1, 1, 0.5,   0,
				  0, 1,   1,   0,
				  1, 2, 1.5, 0.5,
				  0, 1,   1,   1;

	Eigen::VectorXd v_expected(4);
	v_expected << 1, 2, 4, 4;
	
	EXPECT_TRUE(M_expected == M);
	EXPECT_TRUE(v_expected == v);

	//qp.PrintQP_C(std::cout);
}