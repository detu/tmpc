/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

//#include "qp_test_problems.hpp"

#include <tmpc/qp/QpOasesProblem.hpp>
#include <tmpc/qp/Printing.hpp>
#include <tmpc/Matrix.hpp>

#include "gtest_tools_eigen.hpp"

#include <gtest/gtest.h>

#include <array>

using namespace tmpc;

TEST(QpOasesProblemTest, test_MatricesCorrect)
{
	QpOasesProblem p({QpSize(2, 1, 0)});

	StaticMatrix<double, 3, 2, columnMajor> x;

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
	p[0].set_q(DynamicMatrix<double>({{6.}, {7.}}));
	p[0].set_r(DynamicMatrix<double>({{8.}}));

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
	Rand<DynamicMatrix<double>> rand_matrix;
	Rand<DynamicVector<double>> rand_vector;

	for (auto stage = p.begin(); stage != p.end(); ++stage)
	{
		QpSize const sz = stage->size();

		{
			auto tmp = rand_matrix.generate(sz.nx(), sz.nx());
			stage->set_Q(trans(tmp) * tmp);
		}

		{
			auto tmp = rand_matrix.generate(sz.nu(), sz.nu());
			stage->set_R(trans(tmp) * tmp);
		}

		stage->set_S(rand_matrix.generate(sz.nx(), sz.nu()));
		stage->set_q(rand_vector.generate(sz.nx()));
		stage->set_r(rand_vector.generate(sz.nu()));
		stage->set_lbx(rand_vector.generate(sz.nx()));
		stage->set_ubx(rand_vector.generate(sz.nx()));
		stage->set_lbu(rand_vector.generate(sz.nu()));
		stage->set_ubu(rand_vector.generate(sz.nu()));

		if (stage + 1 != p.end())
		{
			stage->set_A(rand_matrix.generate((stage + 1)->size().nx(), sz.nx()));
			stage->set_B(rand_matrix.generate((stage + 1)->size().nx(), sz.nu()));
		}

		stage->set_C(rand_matrix.generate(sz.nc(), sz.nx()));
		stage->set_D(rand_matrix.generate(sz.nc(), sz.nu()));
		stage->set_lbd(rand_vector.generate(sz.nc()));
		stage->set_ubd(rand_vector.generate(sz.nc()));
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

	DynamicMatrix<double> H_expected(total_nx, total_nx, 0.);
	DynamicVector<double> g_expected(total_nx);
	DynamicVector<double> lb_expected(total_nx);
	DynamicVector<double> ub_expected(total_nx);
	DynamicMatrix<double> A_expected(total_nc, total_nx, 0.);
	DynamicVector<double> lbA_expected(total_nc);
	DynamicVector<double> ubA_expected(total_nc);

	size_t i = 0, j = 0, ia = 0;
	for (auto stage = p.begin(); stage != p.end(); ++stage)
	{
		size_t const nx = stage->size().nx();
		size_t const nu = stage->size().nu();
		size_t const nc = stage->size().nc();
		size_t const nx1 = stage + 1 != p.end() ? stage[1].nx() : 0;

		submatrix(H_expected,  i     , j, nx, nx) =       stage->get_Q();	submatrix(H_expected, i     , j + nx, nx, nu) = stage->get_S();
		submatrix(H_expected,  i + nx, j, nu, nx) = trans(stage->get_S());	submatrix(H_expected, i + nx, j + nx, nu, nu) = stage->get_R();

		subvector(g_expected, i, nx) = stage->get_q();	subvector(g_expected, i + nx, nu) = stage->get_r();
		subvector(lb_expected, i, nx) = stage->get_lbx();	subvector(lb_expected, i + nx, nu) = stage->get_lbu();
		subvector(ub_expected, i, nx) = stage->get_ubx();	subvector(ub_expected, i + nx, nu) = stage->get_ubu();

		submatrix(A, ia      , i, nx1, nx) = stage->get_A();	submatrix(A, ia      , i + nx, nx1, nu) = stage->get_B();	submatrix(A, ia, i + nx + nu, nx1, nx1) = -IdentityMatrix<DynamicMatrix<double>>(nx1, 1.);
		submatrix(A, ia + nx1, i, nc , nx) = stage->get_C();	submatrix(A, ia + nx1, i + nx, nx , nu) = stage->get_D();

		subvector(lbA_expected, ia, nx1) = 0.;	subvector(lbA_expected, ia + nx1, nc) = stage->get_lbd();
		subvector(ubA_expected, ia, nx1) = 0.;	subvector(ubA_expected, ia + nx1, nc) = stage->get_ubd();

		i += nx + nu;	j += nx + nu;
		ia += nx + nc;
	}

	EXPECT_EQ(print_wrap(p.H()), print_wrap(H_expected));
	EXPECT_EQ(print_wrap(p.g()), print_wrap(g_expected));
	EXPECT_EQ(print_wrap(p.lb()), print_wrap(lb_expected));
	EXPECT_EQ(print_wrap(p.ub()), print_wrap(ub_expected));
	EXPECT_EQ(print_wrap(p.A()), print_wrap(A_expected));
	EXPECT_EQ(print_wrap(p.lbA()), print_wrap(lbA_expected));
	EXPECT_EQ(print_wrap(p.ubA()), print_wrap(ubA_expected));
}
