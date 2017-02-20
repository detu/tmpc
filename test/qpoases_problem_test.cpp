/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include "gtest_tools_eigen.hpp"
//#include "qp_test_problems.hpp"

#include <tmpc/qp/QpOasesProblem.hpp>
#include <tmpc/qp/Printing.hpp>

#include <gtest/gtest.h>

#include <array>

using namespace tmpc;

TEST(QpOasesProblemTest, test_MatricesCorrect)
{
	tmpc::QpOasesProblem p({QpSize(2, 1, 0)});

	StaticMatrix<double, 3, 2> x;

	// Test H = [Q, S; S', R]
	p[0].set_Q(DynamicMatrix<double>({{1., 2.}, {2., 1.}}));
	p[0].set_R(DynamicMatrix<double>({{3}}));
	p[0].set_S(DynamicMatrix<double>({{4.}, {5.}}));

	EXPECT_EQ(p.H(), DynamicMatrix<double>({
		{1., 2., 4.},
		{2., 1., 5.},
		{4., 5., 3.}
	}));

	// Test g = [q; r]
	p.set_q(0, DynamicMatrix<double>({{6.}, {7.}}));
	p.set_r(0, DynamicMatrix<double>({{8.}}));

	EXPECT_EQ(p.g(), DynamicMatrix<double>({{6.}, {7.}, {8.}}));
}

TEST(QpOasesProblemTest, test_MatricesCorrect_1)
{
	std::array<QpSize, 3> sz = {QpSize(2, 1, 2), QpSize(2, 2, 1), QpSize(3, 1, 3)};
	auto const total_nx = numVariables(sz.begin(), sz.end());
	auto const total_nc = numEqualities(sz.begin(), sz.end()) + numInequalities(sz.begin(), sz.end());

	QpOasesProblem p(sz.begin(), sz.end());

	//---------------
	// Test H, g
	//---------------

	for (auto stage = p.begin(); stage != p.end(); ++stage)
	{
		QpSize const sz = stage->size();

		{
			auto tmp = Rand<DynamicMatrix<double>>(sz.nx(), sz.nx()).generate();
			stage->set_Q(trans(tmp) * tmp);
		}

		{
			auto tmp = Rand<DynamicMatrix<double>>(sz.nu(), sz.nu()).generate();
			stage->set_R(trans(tmp) * tmp);
		}

		stage->set_S(Rand<DynamicMatrix<double>>(sz.nx(), sz.nu()).generate());
		stage->set_q(Rand<DynamicVector<double>>(sz.nx()).generate());
		stage->set_r(Rand<DynamicVector<double>>(sz.nu()).generate());
		stage->set_lbx(Rand<DynamicVector<double>>(sz.nx()).generate());
		stage->set_ubx(Rand<DynamicVector<double>>(sz.nx()).generate());
		stage->set_lbu(Rand<DynamicVector<double>>(sz.nu()).generate());
		stage->set_ubu(Rand<DynamicVector<double>>(sz.nu()).generate());

		if (stage + 1 != p.end())
		{
			stage->set_A(Rand<DynamicMatrix<double>>((stage + 1)->size().nx(), sz.nx()).generate());
			stage->set_B(Rand<DynamicMatrix<double>>((stage + 1)->size().nx(), sz.nu()).generate());
		}

		stage->set_C(Rand<DynamicMatrix<double>>(sz.nc(), sz.nx()).generate());
		stage->set_D(Rand<DynamicMatrix<double>>(sz.nc(), sz.nu()).generate());
		stage->set_lbd(Rand<DynamicVector<double>>(sz.nc()).generate());
		stage->set_ubd(Rand<DynamicVector<double>>(sz.nc()).generate());
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

	// How is the progress in implementing concatenation views in Blaze?
	// https://bitbucket.org/blaze-lib/blaze/issues/44/more-views-concatenation-index-lists

	/*
	EXPECT_EQ(print_wrap(p.H()), print_wrap((MatrixXd <<
		p[0].get_Q(),             p[0].get_S(),              MatrixXd::Zero(sz[0].nx(), sz[1].nx() + sz[1].nu()), MatrixXd::Zero(sz[0].nx(), sz[2].nx() + sz[2].nu()),
		p[0].get_S().transpose(), p[0].get_R(),              MatrixXd::Zero(sz[0].nu(), sz[1].nx() + sz[1].nu()), MatrixXd::Zero(sz[0].nu(), sz[2].nx() + sz[2].nu()),
		MatrixXd::Zero(sz[1].nx(), sz[0].nx() + sz[0].nu()), p[1].get_Q(),             p[1].get_S(),              MatrixXd::Zero(sz[1].nx(), sz[2].nx() + sz[2].nu()),
		MatrixXd::Zero(sz[1].nu(), sz[0].nx() + sz[0].nu()), p[1].get_S().transpose(), p[1].get_R(),              MatrixXd::Zero(sz[1].nu(), sz[2].nx() + sz[2].nu()),
		MatrixXd::Zero(sz[2].nx(), sz[0].nx() + sz[0].nu()), MatrixXd::Zero(sz[2].nx(), sz[1].nx() + sz[1].nu()), p[2].get_Q(),	            p[2].get_S(),
		MatrixXd::Zero(sz[2].nu(), sz[0].nx() + sz[0].nu()), MatrixXd::Zero(sz[2].nu(), sz[1].nx() + sz[1].nu()), p[2].get_S().transpose(),	p[2].get_R()
	).finished()));
	*/
	DynamicMatrix<double> H_expected(total_nx, total_nx, 0.);
	submatrix(H_expected,  0, 0, 2, 2) =       p[0].get_Q();	submatrix(H_expected,  0,  2, 2, 1) = p[0].get_S();
	submatrix(H_expected,  2, 0, 1, 2) = trans(p[0].get_S());	submatrix(H_expected,  2,  2, 1, 1) = p[0].get_R();
	submatrix(H_expected,  3, 3, 2, 2) =       p[1].get_Q();	submatrix(H_expected,  3,  5, 2, 2) = p[1].get_S();
	submatrix(H_expected,  5, 3, 2, 2) = trans(p[1].get_S());	submatrix(H_expected,  5,  5, 2, 2) = p[1].get_R();
	submatrix(H_expected,  7, 7, 3, 3) =       p[2].get_Q();	submatrix(H_expected,  7, 10, 3, 1) = p[2].get_S();
	submatrix(H_expected, 10, 7, 1, 3) = trans(p[2].get_S());	submatrix(H_expected, 10, 10, 1, 1) = p[2].get_R();

	EXPECT_EQ(print_wrap(p.H()), print_wrap(H_expected));

	/*
	EXPECT_EQ(print_wrap(p.g()), print_wrap((VectorXd(total_nx) <<
		p.get_q(0), p.get_r(0), p.get_q(1), p.get_r(1), p.get_q(2), p.get_r(2)
	).finished()));
	*/
	DynamicVector<double> g_expected(total_nx);
	subvector(g_expected, 0, 2) = p[0].get_q();	subvector(g_expected,  2, 1) = p[0].get_r();
	subvector(g_expected, 3, 2) = p[1].get_q();	subvector(g_expected,  5, 2) = p[1].get_r();
	subvector(g_expected, 7, 3) = p[2].get_q();	subvector(g_expected, 10, 1) = p[2].get_r();

	EXPECT_EQ(print_wrap(p.g()), print_wrap(g_expected));


	// *** CONTINUE HERE!!! ***
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
