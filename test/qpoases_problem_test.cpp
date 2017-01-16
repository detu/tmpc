/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include "gtest_tools_eigen.hpp"
//#include "qp_test_problems.hpp"

#include <tmpc/qp/qpOASESProgram.hpp>
#include <tmpc/qp/Printing.hpp>

#include <gtest/gtest.h>

#include <array>

TEST(QpOasesProblemTest, MatricesCorrect_test)
{
	tmpc::qpOASESProgram p({tmpc::QpSize(2, 1, 0)});

	// Test H = [Q, S; S', R]
	p.set_Q(0, (Eigen::MatrixXd(2, 2) << 1., 2., 2., 1).finished());
	p.set_R(0, (Eigen::MatrixXd(1, 1) << 3.).finished());
	p.set_S(0, (Eigen::MatrixXd(2, 1) << 4., 5.).finished());

	EXPECT_PRED2(MatrixApproxEquality(1e-6), p.H(),
		(Eigen::MatrixXd(3, 3) << 1., 2., 4.,
  							      2., 1., 5.,
								  4., 5., 3.).finished());

	// Test g = [q; r]
	p.set_q(0, (Eigen::MatrixXd(2, 1) << 6., 7.).finished());
	p.set_r(0, (Eigen::MatrixXd(1, 1) << 8.).finished());

	EXPECT_PRED2(MatrixApproxEquality(1e-6), p.g(),
			(Eigen::MatrixXd(3, 1) << 6., 7., 8.).finished());
}

TEST(QpOasesProblemTest, H_Problem0_correct)
{
	using tmpc::qpOASESProgram;
	using tmpc::QpSize;
	using Eigen::MatrixXd;
	using Eigen::VectorXd;

	std::array<QpSize, 3> sz = {QpSize(2, 1, 0), QpSize(2, 2, 0), QpSize(3, 1, 0)};
	auto const total_nx = tmpc::numVariables(sz.begin(), sz.end());

	qpOASESProgram p(sz.begin(), sz.end());

	//---------------
	// Test H, g
	//---------------

	for (std::size_t i = 0; i < sz.size(); ++i)
	{
		{
			auto tmp = MatrixXd::Random(sz[i].nx(), sz[i].nx()).eval();
			p.set_Q(i, tmp.transpose() * tmp);
		}

		{
			auto tmp = MatrixXd::Random(sz[i].nu(), sz[i].nu()).eval();
			p.set_R(i, tmp.transpose() * tmp);
		}

		p.set_S(i, MatrixXd::Random(sz[i].nx(), sz[i].nu()));
		p.set_q(i, VectorXd::Random(sz[i].nx()));
		p.set_r(i, VectorXd::Random(sz[i].nu()));
	}

	std::cout << "p.H() = " << std::endl << p.H() << std::endl;
	std::cout << "p.g() = " << std::endl << p.g() << std::endl;

	EXPECT_EQ(print_wrap(p.H()), print_wrap((MatrixXd(total_nx, total_nx) <<
		p.get_Q(0),             p.get_S(0),                  MatrixXd::Zero(sz[0].nx(), sz[1].nx() + sz[1].nu()), MatrixXd::Zero(sz[0].nx(), sz[2].nx() + sz[2].nu()),
		p.get_S(0).transpose(), p.get_R(0),                  MatrixXd::Zero(sz[0].nu(), sz[1].nx() + sz[1].nu()), MatrixXd::Zero(sz[0].nu(), sz[2].nx() + sz[2].nu()),
		MatrixXd::Zero(sz[1].nx(), sz[0].nx() + sz[0].nu()), p.get_Q(1),             p.get_S(1),                  MatrixXd::Zero(sz[1].nx(), sz[2].nx() + sz[2].nu()),
		MatrixXd::Zero(sz[1].nu(), sz[0].nx() + sz[0].nu()), p.get_S(1).transpose(), p.get_R(1),                  MatrixXd::Zero(sz[1].nu(), sz[2].nx() + sz[2].nu()),
		MatrixXd::Zero(sz[2].nx(), sz[0].nx() + sz[0].nu()), MatrixXd::Zero(sz[2].nx(), sz[1].nx() + sz[1].nu()), p.get_Q(2),	            p.get_S(2),
		MatrixXd::Zero(sz[2].nu(), sz[0].nx() + sz[0].nu()), MatrixXd::Zero(sz[2].nu(), sz[1].nx() + sz[1].nu()), p.get_S(2).transpose(),	p.get_R(2)
	).finished()));

	EXPECT_EQ(print_wrap(p.g()), print_wrap((VectorXd(total_nx) <<
		p.get_q(0), p.get_r(0), p.get_q(1), p.get_r(1), p.get_q(2), p.get_r(2)
	).finished()));
}
