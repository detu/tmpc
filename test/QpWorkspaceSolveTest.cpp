/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

//#include <tmpc/qp/CondensingSolver.hpp>
//#include <tmpc/qp/HPMPCSolver.hpp>
#include <tmpc/mpc/MpcQpSize.hpp>	// for mpcQpSize()
#include <tmpc/qp/Printing.hpp>
#include <tmpc/qp/QpOasesWorkspace.hpp>
#include <tmpc/qp/HpmpcWorkspace.hpp>
#include <tmpc/qp/HpipmWorkspace.hpp>
#include <tmpc/Matrix.hpp>

#include "gtest_tools_eigen.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <utility>

using namespace tmpc;

template <typename WS>
class QpWorkspaceSolveTest : public ::testing::Test
{
protected:
	using Workspace = WS;
	using Real = typename Workspace::Real;

	static Workspace problem_0()
	{
		unsigned const NX = 2;
		unsigned const NU = 1;
		unsigned const NZ = NX + NU;
		unsigned const NC = 0;
		unsigned const NCT = 0;
		unsigned const NT = 2;

		typedef StaticMatrix<double, NZ, NZ> StageHessianMatrix;

		const auto sz = mpcQpSize(NT, NX, NU, NC, NCT);
		Workspace qp(sz.begin(), sz.end());
		
		qp[0].lbx(-1.);	qp[0].lbu(-1.);	qp[0].ubx(1.);	qp[0].ubu(1.);
		qp[1].lbx(-1.);	qp[1].lbu(-1.);	qp[1].ubx(1.);	qp[1].ubu(1.);
		qp[2].lbx(-1.);					qp[2].ubx(1.);

		// Stage 0
		StageHessianMatrix H0 {
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9}
		};
		H0 = trans(H0) * H0;	// Make positive definite.

		StaticVector<double, NX> const q0 {0., 0.};
		StaticVector<double, NU> const r0 {0.};

		const DynamicMatrix<double> Q0 = submatrix(H0, 0, 0, NX, NX);
		const DynamicMatrix<double> R0 = submatrix(H0, NX, NX, NU, NU);
		const DynamicMatrix<double> S0 = submatrix(H0, 0, NX, NX, NU);
		const DynamicMatrix<double> S0T = submatrix(H0, NX, 0, NU, NX);

		DynamicMatrix<double> const A0 {{1., 1.}, {0., 1.}};
		DynamicMatrix<double> const B0 {{0.5}, {1.0}};
		DynamicVector<double> a0 {1., 2.};

		// Stage 1
		StageHessianMatrix H1 {
				{1., 2., 3.},
				{4., 5., 6.},
				{7., 8., 9.}
		};
		H1 = trans(H1) * H1;	// Make positive definite.

		StaticVector<double, NX> const q1 {0., 0.};
		StaticVector<double, NU> const r1 {0.};

		const DynamicMatrix<double> Q1 = submatrix(H1, 0, 0, NX, NX);
		const DynamicMatrix<double> R1 = submatrix(H1, NX, NX, NU, NU);
		const DynamicMatrix<double> S1 = submatrix(H1, 0, NX, NX, NU);
		const DynamicMatrix<double> S1T = submatrix(H1, NX, 0, NU, NX);

		DynamicMatrix<double> const A1 {{1., 1.}, {0., 1.}};
		DynamicMatrix<double> const B1 {{0.5}, {1.0}};
		DynamicVector<double> const a1 {1., 2.};

		// Stage 2
		DynamicMatrix<double> H2 {{1., 2.}, {3., 4.}};
		H2 = trans(H2) * H2;	// Make positive definite.

		StaticVector<double, NX> const q2 {0., 0.};

		const DynamicMatrix<double> Q2 = submatrix(H2, 0, 0, NX, NX);

		// Setup QP
		qp[0].Q(Q0);	qp[0].R(R0);	qp[0].S(S0);	qp[0].q(q0);	qp[0].r(r0);
		qp[1].Q(Q1);	qp[1].R(R1);	qp[1].S(S1);	qp[1].q(q1);	qp[1].r(r1);
		qp[2].Q(Q2);									qp[2].q(q2);

		qp[0].A(A0);	
		qp[0].B(B0);		
		qp[0].b(a0);
		
		qp[1].A(A1);	
		qp[1].B(B1);		
		qp[1].b(a1);

		return std::move(qp);
	}

	static Workspace problem_1()
	{
		Workspace qp = problem_0();

		StaticVector<double, 2> x0 {1., 0.};
		qp[0].lbx(x0);	qp[0].ubx(x0);

		qp[0].b(0.);
		qp[1].b(0.);

		return std::move(qp);
	}
};

typedef ::testing::Types<
//		tmpc::CondensingSolver<double>
		QpOasesWorkspace,
		HpmpcWorkspace<double>,
		HpipmWorkspace<double>
	> SolverTypes;

TYPED_TEST_CASE(QpWorkspaceSolveTest, SolverTypes);

/// \brief Check if QPSolver move constructor works and the solver works after move constructor.
TYPED_TEST(QpWorkspaceSolveTest, testMoveConstructor)
{
	auto ws = TestFixture::problem_0();
	auto ws1 = std::move(ws);

	ws1.solve();

	using Vector = DynamicVector<typename TestFixture::Real>;

	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws1[0].x(), (Vector {1., -1.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws1[0].u(), (Vector {-1.}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws1[1].x(), (Vector {0.5, 0.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws1[1].u(), (Vector {-1.}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws1[2].x(), (Vector {1., 1.}));
}

TYPED_TEST(QpWorkspaceSolveTest, testSolve0)
{
	auto ws = TestFixture::problem_0();

	ws.solve();

	using Vector = DynamicVector<typename TestFixture::Real>;

	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[0].x(), (Vector {1., -1.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[0].u(), (Vector {-1.}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[1].x(), (Vector {0.5, 0.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[1].u(), (Vector {-1.}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[2].x(), (Vector {1., 1.}));
}

TYPED_TEST(QpWorkspaceSolveTest, testSolve1)
{
	auto ws = TestFixture::problem_1();

	ws.solve();

	using Vector = DynamicVector<typename TestFixture::Real>;

	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[0].x(), (Vector {1., 0.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[0].u(), (Vector {-0.690877362606266}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[1].x(), (Vector {0.654561318696867, -0.690877362606266}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[1].u(), (Vector {0.215679569867116}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), ws[2].x(), (Vector {0.0715237410241597, -0.475197792739149}));
}
