#pragma once

#include <tmpc/mpc/MpcQpSize.hpp>	// for mpcQpSize()
#include <tmpc/qp/Printing.hpp>

#include <tmpc/Matrix.hpp>

#include <tmpc/qp/Printing.hpp>

#include "../gtest_tools_eigen.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <fstream>

namespace tmpc :: testing
{
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
			Workspace ws(sz.begin(), sz.end());
			auto qp = ws.problem();
			
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

			return std::move(ws);
		}

		static Workspace problem_1()
		{
			Workspace ws = problem_0();
			auto qp = ws.problem();

			StaticVector<double, 2> x0 {1., 0.};
			qp[0].lbx(x0);	qp[0].ubx(x0);

			qp[0].b(0.);
			qp[1].b(0.);

			return std::move(ws);
		}
	};

	TYPED_TEST_CASE_P(QpWorkspaceSolveTest);

	/// \brief Check if QPSolver move constructor works and the solver works after move constructor.
	TYPED_TEST_P(QpWorkspaceSolveTest, testMoveConstructor)
	{
		auto ws = TestFixture::problem_0();
		auto ws1 = std::move(ws);

		ws1.solve();
		auto const sol = ws1.solution();

		using Vector = DynamicVector<typename TestFixture::Real>;

		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[0].x(), (Vector {1., -1.}));
		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[0].u(), (Vector {-1.}));

		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[1].x(), (Vector {0.5, 0.}));
		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[1].u(), (Vector {-1.}));

		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[2].x(), (Vector {1., 1.}));
	}

	TYPED_TEST_P(QpWorkspaceSolveTest, testSolve0)
	{
		auto ws = TestFixture::problem_0();

		ws.solve();
		auto const sol = ws.solution();

		using Vector = DynamicVector<typename TestFixture::Real>;

		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[0].x(), (Vector {1., -1.}));
		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[0].u(), (Vector {-1.}));

		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[1].x(), (Vector {0.5, 0.}));
		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[1].u(), (Vector {-1.}));

		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[2].x(), (Vector {1., 1.}));
	}

	TYPED_TEST_P(QpWorkspaceSolveTest, testSolve1)
	{
		auto ws = TestFixture::problem_1();

		ws.solve();
		auto const sol = ws.solution();

		using Vector = DynamicVector<typename TestFixture::Real>;

		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[0].x(), (Vector {1., 0.}));
		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[0].u(), (Vector {-0.690877362606266}));

		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[1].x(), (Vector {0.654561318696867, -0.690877362606266}));
		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[1].u(), (Vector {0.215679569867116}));

		EXPECT_PRED2(MatrixApproxEquality(1e-6), sol[2].x(), (Vector {0.0715237410241597, -0.475197792739149}));
	}

	TYPED_TEST_P(QpWorkspaceSolveTest, DISABLED_testSolve1stage1d)
	{
		typename TestFixture::Workspace ws { std::array<OpSize, 1> { OpSize(1, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(2.);
		problem[0].q(-1.);
		problem[0].lbx(-100.);
		problem[0].ubx(100.);

		ws.solve();
		auto solution = ws.solution();

		using Vector = DynamicVector<typename TestFixture::Real>;

		EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].x(), (Vector {0.5}));
	}

	///
	/// minimize (1/2 * Q0 * x0^2 - q0 * x0) + (1/2 * Q1 * x1^2 - q1 * x1)
	/// s.t. x1 = A0 * x0 + b0
	///      -100 <= x0 <= 100
	///      -100 <= x1 <= 100
	///
	TYPED_TEST_P(QpWorkspaceSolveTest, testSolve2stage1d)
	{
		using Real = typename TestFixture::Real;

		Real const Q0 = 2., q0 = -1., Q1 = 1., q1 = 0., A0 = 0.2, b0 = 1.;

		typename TestFixture::Workspace ws {std::array<OpSize, 2> { OpSize(1, 0, 0), OpSize(1, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(Q0);
		problem[0].q(q0);
		problem[0].A(A0);
		problem[0].b(b0);
		problem[0].lbx(-100.);
		problem[0].ubx(100.);

		problem[1].Q(Q1);
		problem[1].q(q1);
		problem[1].lbx(-100.);
		problem[1].ubx(100.);

		ws.solve();

		/*
		try
		{
			ws.solve();
		}
		catch (HpmpcException const&)
		{
			auto const& hpmpc_ws = reinterpret_cast<HpmpcWorkspace<double> const&>(ws);
			for (auto const& stat : hpmpc_ws.stat())
			{
				for (auto val : stat)
					std::cout << val << "\t";
				std::cout << std::endl;
			}

			{
				std::ofstream failed_qp("failed_qp.m");
				PrintMultistageQpMatlab(failed_qp, ws.problem(), "qp");
			}

			throw;
		}
		*/

		auto solution = ws.solution();

		using Vector = DynamicVector<Real>;

		Real const x0_opt = -(q0 + Q1 * A0 * b0 + q1 * A0) / (Q0 + Q1 * A0 * A0);

		EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].x(), (Vector {x0_opt}));
		EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[1].x(), (Vector {A0 * x0_opt + b0}));
	}

	/**
	 \brief This QP has 3 states and 2 inputs, so that the S matrix is 3x2.
	 
	 It checks if the element order of S (row-major or col-major) is correct 
	 for the solvers that depend on it (like hpmpc or hpipm).
	 */
	TYPED_TEST_P(QpWorkspaceSolveTest, testSolve2stage5d)
	{
		using Real = typename TestFixture::Real;

		typename TestFixture::Workspace ws {std::array<OpSize, 2> { OpSize(3, 2, 0), OpSize(0, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(DynamicMatrix<Real>({
			{2.73449304081811,	1.88589647487061,	2.07850896302612},
			{1.88589647487061,	2.23398400982993,	2.04607087035366},
			{2.07850896302612,	2.04607087035366,	2.75909914758841}
		}));
		
		problem[0].R(DynamicMatrix<Real>({
			{2.58480403401811,	2.27683785836085},
			{2.27683785836085,	2.48531656389865}
		}));

		problem[0].S(DynamicMatrix<Real>({
			{1.94423942432824,	1.95671439063179},
			{2.31644542710892,	2.08748866537729},
			{2.46061457306194,	1.94728938270016}
		}));

		problem[0].q(DynamicVector<Real>({
			0.276025076998578,
			0.679702676853675,
			0.655098003973841
		}));

		problem[0].r(DynamicVector<Real>({
			0.162611735194631,
			0.118997681558377
		}));

		problem[0].A(StaticMatrix<Real, 0, 3>(0.));
		problem[0].B(StaticMatrix<Real, 0, 2>(0.));
		problem[0].b(StaticVector<Real, 0>(0.));
		problem[0].lbx(-10000.);
		problem[0].ubx(10000.);
		problem[0].lbu(-10000.);
		problem[0].ubu(10000.);

		try
		{
			ws.solve();
			auto solution = ws.solution();

			using Vector = DynamicVector<Real>;

			EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].x(), (Vector {
				146.566682434017,
				-427.218558989821,
				-345.347969700289
			}));

			EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].u(), (Vector {
				769.663140469139,
				-191.122524763114
			}));

			EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[1].x(), (Vector {}));
			EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[1].u(), (Vector {}));
		}
		catch (std::runtime_error const&)
		{
			std::cout << "------ PROBLEM ------" << std::endl;
			for (auto const& p : ws.problem())
				std::cout << p << std::endl;

			std::cout << "------ SOLUTION ------" << std::endl;
			for (auto const& s : ws.solution())
				std::cout << s << std::endl;

			throw;
		}
	}

	REGISTER_TYPED_TEST_CASE_P(QpWorkspaceSolveTest,
        testMoveConstructor, testSolve0, testSolve1, DISABLED_testSolve1stage1d, testSolve2stage1d, testSolve2stage5d);
}