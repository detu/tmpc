/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include <tmpc/test_tools.hpp>
//#include "qp_test_problems.hpp"

#include <tmpc/qp/QpOasesProblem.hpp>
#include <tmpc/qp/Printing.hpp>

#include <gtest/gtest.h>

#include <array>

TEST(QpOasesProblemTest, test_MatricesCorrect)
{
	tmpc::QpOasesProblem p({tmpc::QpSize(2, 1, 0)});

	// Test H = [Q, S; S', R]
	p.set_Q(0, (Eigen::MatrixXd(2, 2) << 1., 2., 2., 1).finished());
	p.set_R(0, (Eigen::MatrixXd(1, 1) << 3.).finished());
	p.set_S(0, (Eigen::MatrixXd(2, 1) << 4., 5.).finished());

	EXPECT_EQ(p.H(),
		(Eigen::MatrixXd(3, 3) << 1., 2., 4.,
  							      2., 1., 5.,
								  4., 5., 3.).finished());

	// Test g = [q; r]
	p.set_q(0, (Eigen::MatrixXd(2, 1) << 6., 7.).finished());
	p.set_r(0, (Eigen::MatrixXd(1, 1) << 8.).finished());

	EXPECT_EQ(p.g(), (Eigen::MatrixXd(3, 1) << 6., 7., 8.).finished());
}

TEST(QpOasesProblemTest, test_MatricesCorrect_1)
{
	using tmpc::QpOasesProblem;
	using tmpc::QpSize;
	using Eigen::MatrixXd;
	using Eigen::VectorXd;

	std::array<QpSize, 3> sz = {QpSize(2, 1, 2), QpSize(2, 2, 1), QpSize(3, 1, 3)};
	auto const total_nx = tmpc::numVariables(sz.begin(), sz.end());
	auto const total_nc = tmpc::numEqualities(sz.begin(), sz.end()) + tmpc::numInequalities(sz.begin(), sz.end());

	QpOasesProblem p(sz.begin(), sz.end());

	//---------------
	// Test H, g
	//---------------

	for (auto stage = p.begin(); stage != p.end(); ++stage)
	{
		QpSize const sz = stage->size();

		{
			auto tmp = MatrixXd::Random(sz.nx(), sz.nx()).eval();
			stage->set_Q(tmp.transpose() * tmp);
		}

		{
			auto tmp = MatrixXd::Random(sz.nu(), sz.nu()).eval();
			stage->set_R(tmp.transpose() * tmp);
		}

		stage->set_S(MatrixXd::Random(sz.nx(), sz.nu()));
		stage->set_q(VectorXd::Random(sz.nx()));
		stage->set_r(VectorXd::Random(sz.nu()));
		stage->set_lbx(VectorXd::Random(sz.nx()));
		stage->set_ubx(VectorXd::Random(sz.nx()));
		stage->set_lbu(VectorXd::Random(sz.nu()));
		stage->set_ubu(VectorXd::Random(sz.nu()));

		if (stage + 1 != p.end())
		{
			stage->set_A(MatrixXd::Random((stage + 1)->size().nx(), sz.nx()));
			stage->set_B(MatrixXd::Random((stage + 1)->size().nx(), sz.nu()));
		}

		stage->set_C(MatrixXd::Random(sz.nc(), sz.nx()));
		stage->set_D(MatrixXd::Random(sz.nc(), sz.nu()));
		stage->set_lbd(VectorXd::Random(sz.nc()));
		stage->set_ubd(VectorXd::Random(sz.nc()));
	}

	/*
	std::cout << "p.H() = " << std::endl << p.H() << std::endl;
	std::cout << "p.g() = " << std::endl << p.g() << std::endl;
	std::cout << "p.lb() = " << std::endl << p.lb() << std::endl;
	std::cout << "p.ub() = " << std::endl << p.ub() << std::endl;
	std::cout << "p.A() = " << std::endl << p.A() << std::endl;
	std::cout << "p.lbA() = " << std::endl << p.lbA() << std::endl;
	std::cout << "p.ubA() = " << std::endl << p.ubA() << std::endl;
	*/

	EXPECT_EQ(print_wrap(p.H()), print_wrap((MatrixXd(total_nx, total_nx) <<
		p[0].get_Q(),             p[0].get_S(),              MatrixXd::Zero(sz[0].nx(), sz[1].nx() + sz[1].nu()), MatrixXd::Zero(sz[0].nx(), sz[2].nx() + sz[2].nu()),
		p[0].get_S().transpose(), p[0].get_R(),              MatrixXd::Zero(sz[0].nu(), sz[1].nx() + sz[1].nu()), MatrixXd::Zero(sz[0].nu(), sz[2].nx() + sz[2].nu()),
		MatrixXd::Zero(sz[1].nx(), sz[0].nx() + sz[0].nu()), p[1].get_Q(),             p[1].get_S(),              MatrixXd::Zero(sz[1].nx(), sz[2].nx() + sz[2].nu()),
		MatrixXd::Zero(sz[1].nu(), sz[0].nx() + sz[0].nu()), p[1].get_S().transpose(), p[1].get_R(),              MatrixXd::Zero(sz[1].nu(), sz[2].nx() + sz[2].nu()),
		MatrixXd::Zero(sz[2].nx(), sz[0].nx() + sz[0].nu()), MatrixXd::Zero(sz[2].nx(), sz[1].nx() + sz[1].nu()), p[2].get_Q(),	            p[2].get_S(),
		MatrixXd::Zero(sz[2].nu(), sz[0].nx() + sz[0].nu()), MatrixXd::Zero(sz[2].nu(), sz[1].nx() + sz[1].nu()), p[2].get_S().transpose(),	p[2].get_R()
	).finished()));

	EXPECT_EQ(print_wrap(p.g()), print_wrap((VectorXd(total_nx) <<
		p.get_q(0), p.get_r(0), p.get_q(1), p.get_r(1), p.get_q(2), p.get_r(2)
	).finished()));

	EXPECT_EQ(print_wrap(p.lb()), print_wrap((VectorXd(total_nx) <<
		p[0].get_lbx(), p[0].get_lbu(), p[1].get_lbx(), p[1].get_lbu(), p[2].get_lbx(), p[2].get_lbu()
	).finished()));

	EXPECT_EQ(print_wrap(p.ub()), print_wrap((VectorXd(total_nx) <<
		p[0].get_ubx(), p[0].get_ubu(), p[1].get_ubx(), p[1].get_ubu(), p[2].get_ubx(), p[2].get_ubu()
	).finished()));

	EXPECT_EQ(print_wrap(p.A()), print_wrap((MatrixXd(total_nc, total_nx) <<
		p[0].get_A(),                           p[0].get_B(),                           -MatrixXd::Identity(sz[1].nx(), sz[1].nx()), MatrixXd::Zero(sz[1].nx(), sz[1].nu()),  MatrixXd::Zero(sz[1].nx(), sz[2].nx()),     MatrixXd::Zero(sz[1].nx(), sz[2].nu()),
		p[0].get_C(),                           p[0].get_D(),                            MatrixXd::Zero(sz[0].nc(), sz[1].nx()),     MatrixXd::Zero(sz[0].nc(), sz[1].nu()),  MatrixXd::Zero(sz[0].nc(), sz[2].nx()),     MatrixXd::Zero(sz[0].nc(), sz[2].nu()),
		MatrixXd::Zero(sz[2].nx(), sz[0].nx()), MatrixXd::Zero(sz[2].nx(), sz[0].nu()),  p[1].get_A(),                               p[1].get_B(),                           -MatrixXd::Identity(sz[2].nx(), sz[2].nx()), MatrixXd::Zero(sz[2].nx(), sz[2].nu()),
	    MatrixXd::Zero(sz[1].nc(), sz[0].nx()), MatrixXd::Zero(sz[1].nc(), sz[0].nu()),  p[1].get_C(),                               p[1].get_D(),                            MatrixXd::Zero(sz[1].nc(), sz[2].nx()),     MatrixXd::Zero(sz[1].nc(), sz[2].nu()),
		MatrixXd::Zero(sz[2].nc(), sz[0].nx()), MatrixXd::Zero(sz[2].nc(), sz[0].nu()),  MatrixXd::Zero(sz[2].nc(), sz[1].nx()),     MatrixXd::Zero(sz[2].nc(), sz[1].nu()),  p[2].get_C(),                               p[2].get_D()
	).finished()));

	EXPECT_EQ(print_wrap(p.lbA()), print_wrap((VectorXd(total_nc) <<
		VectorXd::Zero(sz[1].nx()), p[0].get_lbd(), VectorXd::Zero(sz[2].nx()), p[1].get_lbd(), p[2].get_lbd()
	).finished()));

	EXPECT_EQ(print_wrap(p.ubA()), print_wrap((VectorXd(total_nc) <<
		VectorXd::Zero(sz[1].nx()), p[0].get_ubd(), VectorXd::Zero(sz[2].nx()), p[1].get_ubd(), p[2].get_ubd()
	).finished()));
}
