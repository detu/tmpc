/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

//#include "qp_test_problems.hpp"

#include <tmpc/qp/QpOasesWorkspace.hpp>
#include <tmpc/qp/Printing.hpp>
#include <tmpc/Matrix.hpp>

#include "gtest_tools_eigen.hpp"

#include <gtest/gtest.h>

#include <array>

using namespace tmpc;

TEST(QpOasesWorkspaceTest, testMatricesCorrect)
{
	QpOasesWorkspace ws(std::array<QpSize, 1> {QpSize{2, 1, 0}});
	auto p = ws.problem();

	StaticMatrix<double, 3, 2, columnMajor> x;

	// Test H = [Q, S; S', R]
	p[0].Q(DynamicMatrix<double>({{1., 2.}, {2., 1.}}));
	p[0].R(DynamicMatrix<double>({{3}}));
	p[0].S(DynamicMatrix<double>({{4.}, {5.}}));

	EXPECT_EQ(ws.H(), DynamicMatrix<double>({
		{1., 2., 4.},
		{2., 1., 5.},
		{4., 5., 3.}
	}));

	// Test g = [q; r]
	p[0].q(DynamicMatrix<double>({{6.}, {7.}}));
	p[0].r(DynamicMatrix<double>({{8.}}));

	EXPECT_EQ(ws.g(), DynamicMatrix<double>({{6.}, {7.}, {8.}}));
}

TEST(QpOasesWorkspaceTest, testMatricesCorrect1)
{
	std::array<QpSize, 3> sz = {QpSize(2, 1, 2), QpSize(2, 2, 1), QpSize(3, 1, 3)};
	auto const total_nx = numVariables(sz.begin(), sz.end());
	auto const total_nc = numEqualities(sz.begin(), sz.end()) + numInequalities(sz.begin(), sz.end());

	QpOasesWorkspace ws(sz.begin(), sz.end());
	auto p = ws.problem();

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
			stage->Q(trans(tmp) * tmp);
		}

		{
			auto tmp = rand_matrix.generate(sz.nu(), sz.nu());
			stage->R(trans(tmp) * tmp);
		}

		stage->S(rand_matrix.generate(sz.nx(), sz.nu()));
		stage->q(rand_vector.generate(sz.nx()));
		stage->r(rand_vector.generate(sz.nu()));
		stage->lbx(rand_vector.generate(sz.nx()));
		stage->ubx(rand_vector.generate(sz.nx()));
		stage->lbu(rand_vector.generate(sz.nu()));
		stage->ubu(rand_vector.generate(sz.nu()));

		if (stage + 1 != p.end())
		{
			stage->A(rand_matrix.generate((stage + 1)->size().nx(), sz.nx()));
			stage->B(rand_matrix.generate((stage + 1)->size().nx(), sz.nu()));
		}

		stage->C(rand_matrix.generate(sz.nc(), sz.nx()));
		stage->D(rand_matrix.generate(sz.nc(), sz.nu()));
		stage->lbd(rand_vector.generate(sz.nc()));
		stage->ubd(rand_vector.generate(sz.nc()));
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

	auto constexpr nan = std::numeric_limits<double>::signaling_NaN();
	DynamicMatrix<double> H_expected(total_nx, total_nx, 0.);
	DynamicVector<double> g_expected(total_nx, nan);
	DynamicVector<double> lb_expected(total_nx, nan);
	DynamicVector<double> ub_expected(total_nx, nan);
	DynamicMatrix<double> A_expected(total_nc, total_nx, 0.);
	DynamicVector<double> lbA_expected(total_nc, nan);
	DynamicVector<double> ubA_expected(total_nc, nan);

	size_t i = 0, ia = 0;
	for (auto stage = p.begin(); stage != p.end(); ++stage)
	{
		size_t const nx = stage->size().nx();
		size_t const nu = stage->size().nu();
		size_t const nc = stage->size().nc();
		size_t const nx1 = stage + 1 != p.end() ? stage[1].size().nx() : 0;

		submatrix(H_expected,  i     , i, nx, nx) =       stage->Q();	submatrix(H_expected, i     , i + nx, nx, nu) = stage->S();
		submatrix(H_expected,  i + nx, i, nu, nx) = trans(stage->S());	submatrix(H_expected, i + nx, i + nx, nu, nu) = stage->R();

		subvector(g_expected, i, nx) = stage->q();	subvector(g_expected, i + nx, nu) = stage->r();
		subvector(lb_expected, i, nx) = stage->lbx();	subvector(lb_expected, i + nx, nu) = stage->lbu();
		subvector(ub_expected, i, nx) = stage->ubx();	subvector(ub_expected, i + nx, nu) = stage->ubu();

		submatrix(A_expected, ia      , i, nx1, nx) = stage->A();	submatrix(A_expected, ia      , i + nx, nx1, nu) = stage->B();	submatrix(A_expected, ia, i + nx + nu, nx1, nx1) = -IdentityMatrix<double>(nx1);
		submatrix(A_expected, ia + nx1, i, nc , nx) = stage->C();	submatrix(A_expected, ia + nx1, i + nx, nc , nu) = stage->D();

		subvector(lbA_expected, ia, nx1) = 0.;	subvector(lbA_expected, ia + nx1, nc) = stage->lbd();
		subvector(ubA_expected, ia, nx1) = 0.;	subvector(ubA_expected, ia + nx1, nc) = stage->ubd();

		i += nx + nu;
		ia += nx1 + nc;
	}

	ASSERT_EQ(ia, rows(A_expected));

	EXPECT_EQ(print_wrap(ws.H()), print_wrap(H_expected));
	EXPECT_EQ(print_wrap(ws.g()), print_wrap(g_expected));
	EXPECT_EQ(print_wrap(ws.lb()), print_wrap(lb_expected));
	EXPECT_EQ(print_wrap(ws.ub()), print_wrap(ub_expected));
	EXPECT_EQ(print_wrap(ws.A()), print_wrap(A_expected));
	EXPECT_EQ(print_wrap(ws.lbA()), print_wrap(lbA_expected));
	EXPECT_EQ(print_wrap(ws.ubA()), print_wrap(ubA_expected));
}
