#pragma once

#include <tmpc/mpc/MpcOcpSize.hpp>	// for mpcOcpSize()

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/Testing.hpp>

#include <iostream>
#include <fstream>


namespace tmpc :: testing
{
	template <typename WS>
	class QpWorkspaceSolveTest 
	: 	public Test
	{
	protected:
		using Workspace = WS;
		using Kernel = typename WS::Kernel;
		using Real = typename Workspace::Real;
		using Vector = DynamicVector<Kernel>;
		using Matrix = DynamicMatrix<Kernel>;

		static Workspace problem_0()
		{
			unsigned const NX = 2;
			unsigned const NU = 1;
			unsigned const NZ = NX + NU;
			unsigned const NC = 0;
			unsigned const NCT = 0;
			unsigned const NT = 2;

			typedef StaticMatrix<Kernel, NZ, NZ> StageHessianMatrix;

			const auto sz = mpcOcpSize(NT, NX, NU, NC, NCT);
			Workspace ws(sz.begin(), sz.end());
			auto qp = ws.problem();
			
			qp[0].lbx(-1.);	qp[0].lbu(-1.);	qp[0].ubx(1.);	qp[0].ubu(1.);
			qp[1].lbx(-1.);	qp[1].lbu(-1.);	qp[1].ubx(1.);	qp[1].ubu(1.);
			qp[2].lbx(-1.);					qp[2].ubx(1.);

			// Stage 0
			StageHessianMatrix H0 {
				{67,   78,   90},
				{78,   94,  108},
				{90,  108,  127}
			};

			StaticVector<Kernel, NX> const q0 {0., 0.};
			StaticVector<Kernel, NU> const r0 {0.};

			const DynamicMatrix<Kernel> Q0 = submatrix(H0, 0, 0, NX, NX);
			const DynamicMatrix<Kernel> R0 = submatrix(H0, NX, NX, NU, NU);
			const DynamicMatrix<Kernel> S0 = submatrix(H0, 0, NX, NX, NU);
			const DynamicMatrix<Kernel> S0T = submatrix(H0, NX, 0, NU, NX);

			DynamicMatrix<Kernel> const A0 {{1., 1.}, {0., 1.}};
			DynamicMatrix<Kernel> const B0 {{0.5}, {1.0}};
			Vector a0 {1., 2.};

			// Stage 1
			StageHessianMatrix H1 {
				{67,   78,   90},
				{78,   94,  108},
				{90,  108,  127}
			};

			StaticVector<Kernel, NX> const q1 {0., 0.};
			StaticVector<Kernel, NU> const r1 {0.};

			const DynamicMatrix<Kernel> Q1 = submatrix(H1, 0, 0, NX, NX);
			const DynamicMatrix<Kernel> R1 = submatrix(H1, NX, NX, NU, NU);
			const DynamicMatrix<Kernel> S1 = submatrix(H1, 0, NX, NX, NU);
			const DynamicMatrix<Kernel> S1T = submatrix(H1, NX, 0, NU, NX);

			DynamicMatrix<Kernel> const A1 {{1., 1.}, {0., 1.}};
			DynamicMatrix<Kernel> const B1 {{0.5}, {1.0}};
			Vector const a1 {1., 2.};

			// Stage 2
			DynamicMatrix<Kernel> H2 {{1., 2.}, {3., 4.}};
			H2 = trans(H2) * H2;	// Make positive definite.

			StaticVector<Kernel, NX> const q2 {0., 0.};

			const DynamicMatrix<Kernel> Q2 = submatrix(H2, 0, 0, NX, NX);

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

			StaticVector<Kernel, 2> x0 {1., 0.};
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

		EXPECT_TRUE(approxEqual(sol[0].x(), (DynamicVector<typename TestFixture::Kernel> {1., -1.}), 1e-6));
		EXPECT_TRUE(approxEqual(sol[0].u(), (DynamicVector<typename TestFixture::Kernel> {-1.}), 1e-6));

		EXPECT_TRUE(approxEqual(sol[1].x(), (DynamicVector<typename TestFixture::Kernel> {0.5, 0.}), 1e-6));
		EXPECT_TRUE(approxEqual(sol[1].u(), (DynamicVector<typename TestFixture::Kernel> {-1.}), 1e-6));

		EXPECT_TRUE(approxEqual(sol[2].x(), (DynamicVector<typename TestFixture::Kernel> {1., 1.}), 1e-6));
	}

	TYPED_TEST_P(QpWorkspaceSolveTest, testSolve0)
	{
		auto ws = TestFixture::problem_0();

		ws.solve();
		auto const sol = ws.solution();

		EXPECT_TRUE(approxEqual(sol[0].x(), (DynamicVector<typename TestFixture::Kernel> {1., -1.}), 1e-6));
		EXPECT_TRUE(approxEqual(sol[0].u(), (DynamicVector<typename TestFixture::Kernel> {-1.}), 1e-6));

		EXPECT_TRUE(approxEqual(sol[1].x(), (DynamicVector<typename TestFixture::Kernel> {0.5, 0.}), 1e-6));
		EXPECT_TRUE(approxEqual(sol[1].u(), (DynamicVector<typename TestFixture::Kernel> {-1.}), 1e-6));

		EXPECT_TRUE(approxEqual(sol[2].x(), (DynamicVector<typename TestFixture::Kernel> {1., 1.}), 1e-6));
	}

	TYPED_TEST_P(QpWorkspaceSolveTest, testSolve1)
	{
		auto ws = TestFixture::problem_1();

		ws.solve();
		auto const sol = ws.solution();

		EXPECT_TRUE(approxEqual(sol[0].x(), (DynamicVector<typename TestFixture::Kernel> {1., 0.}), 1e-6));
		EXPECT_TRUE(approxEqual(sol[0].u(), (DynamicVector<typename TestFixture::Kernel> {-0.68098253759615734}), 1e-6));

		EXPECT_TRUE(approxEqual(sol[1].x(), (DynamicVector<typename TestFixture::Kernel> {0.65950873120185627, -0.68098253759609839}), 1e-6));
		EXPECT_TRUE(approxEqual(sol[1].u(), (DynamicVector<typename TestFixture::Kernel> {0.20174225742383531}), 1e-6));

		EXPECT_TRUE(approxEqual(sol[2].x(), (DynamicVector<typename TestFixture::Kernel> {0.079397322317675517, -0.47924028017226311}), 1e-6));
	}
	

	TYPED_TEST_P(QpWorkspaceSolveTest, DISABLED_testSolve1stage1d)
	{
		typename TestFixture::Workspace ws { std::array<OcpSize, 1> { OcpSize(1, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(2.);
		problem[0].q(-1.);
		problem[0].lbx(-100.);
		problem[0].ubx(100.);

		ws.solve();
		auto solution = ws.solution();

		EXPECT_TRUE(approxEqual(solution[0].x(), (DynamicVector<typename TestFixture::Kernel> {0.5}), 1e-6));
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

		typename TestFixture::Workspace ws {std::array<OcpSize, 2> { OcpSize(1, 0, 0), OcpSize(1, 0, 0) } };

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
		auto solution = ws.solution();

		Real const x0_opt = -(q0 + Q1 * A0 * b0 + q1 * A0) / (Q0 + Q1 * A0 * A0);

		EXPECT_TRUE(approxEqual(solution[0].x(), (DynamicVector<typename TestFixture::Kernel> {x0_opt}), 1e-6));
		EXPECT_TRUE(approxEqual(solution[1].x(), (DynamicVector<typename TestFixture::Kernel> {A0 * x0_opt + b0}), 1e-6));
	}

	/**
	 \brief This QP has 3 states and 2 inputs, so that the S matrix is 3x2.
	 
	 It checks if the element order of S (row-major or col-major) is correct 
	 for the solvers that depend on it (like hpmpc or hpipm).
	 */
	TYPED_TEST_P(QpWorkspaceSolveTest, testSolve2stage5d)
	{
		using Real = typename TestFixture::Real;

		typename TestFixture::Workspace ws {std::array<OcpSize, 2> { OcpSize(3, 2, 0), OcpSize(0, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(DynamicMatrix<typename TestFixture::Kernel>({
			{2.73449304081811,	1.88589647487061,	2.07850896302612},
			{1.88589647487061,	2.23398400982993,	2.04607087035366},
			{2.07850896302612,	2.04607087035366,	2.75909914758841}
		}));
		
		problem[0].R(DynamicMatrix<typename TestFixture::Kernel>({
			{2.58480403401811,	2.27683785836085},
			{2.27683785836085,	2.48531656389865}
		}));

		problem[0].S(DynamicMatrix<typename TestFixture::Kernel>({
			{1.94423942432824,	1.95671439063179},
			{2.31644542710892,	2.08748866537729},
			{2.46061457306194,	1.94728938270016}
		}));

		problem[0].q(DynamicVector<typename TestFixture::Kernel>({
			0.276025076998578,
			0.679702676853675,
			0.655098003973841
		}));

		problem[0].r(DynamicVector<typename TestFixture::Kernel>({
			0.162611735194631,
			0.118997681558377
		}));

		problem[0].A(StaticMatrix<typename TestFixture::Kernel, 0, 3>(0.));
		problem[0].B(StaticMatrix<typename TestFixture::Kernel, 0, 2>(0.));
		problem[0].b(StaticVector<typename TestFixture::Kernel, 0>(0.));
		problem[0].lbx(-10000.);
		problem[0].ubx(10000.);
		problem[0].lbu(-10000.);
		problem[0].ubu(10000.);

		ws.solve();
		auto solution = ws.solution();

		EXPECT_TRUE(approxEqual(solution[0].x(), (DynamicVector<typename TestFixture::Kernel> {
			146.566682434017,
			-427.218558989821,
			-345.347969700289
		}), 1e-6));

		EXPECT_TRUE(approxEqual(solution[0].u(), (DynamicVector<typename TestFixture::Kernel> {
			769.663140469139,
			-191.122524763114
		}), 1e-6));

		EXPECT_TRUE(approxEqual(solution[1].x(), (DynamicVector<typename TestFixture::Kernel> {}), 1e-6));
		EXPECT_TRUE(approxEqual(solution[1].u(), (DynamicVector<typename TestFixture::Kernel> {}), 1e-6));
	}

	/**
	 \brief Check that the case when all state and input bounds are infinite is correctly handled by a QP solver.
	 */
	TYPED_TEST_P(QpWorkspaceSolveTest, testInfiniteBoundsAll)
	{
		using Real = typename TestFixture::Real;

		typename TestFixture::Workspace ws {std::array<OcpSize, 2> { OcpSize(3, 2, 0), OcpSize(0, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(DynamicMatrix<typename TestFixture::Kernel>({
			{2.73449304081811,	1.88589647487061,	2.07850896302612},
			{1.88589647487061,	2.23398400982993,	2.04607087035366},
			{2.07850896302612,	2.04607087035366,	2.75909914758841}
		}));
		
		problem[0].R(DynamicMatrix<typename TestFixture::Kernel>({
			{2.58480403401811,	2.27683785836085},
			{2.27683785836085,	2.48531656389865}
		}));

		problem[0].S(DynamicMatrix<typename TestFixture::Kernel>({
			{1.94423942432824,	1.95671439063179},
			{2.31644542710892,	2.08748866537729},
			{2.46061457306194,	1.94728938270016}
		}));

		problem[0].q(DynamicVector<typename TestFixture::Kernel>({
			0.276025076998578,
			0.679702676853675,
			0.655098003973841
		}));

		problem[0].r(DynamicVector<typename TestFixture::Kernel>({
			0.162611735194631,
			0.118997681558377
		}));

		problem[0].A(StaticMatrix<typename TestFixture::Kernel, 0, 3>(0.));
		problem[0].B(StaticMatrix<typename TestFixture::Kernel, 0, 2>(0.));
		problem[0].b(StaticVector<typename TestFixture::Kernel, 0>(0.));
		problem[0].lbx(-inf<Real>());
		problem[0].ubx(inf<Real>());
		problem[0].lbu(-inf<Real>());
		problem[0].ubu(inf<Real>());

		ws.solve();
		auto solution = ws.solution();

		using Vector = DynamicVector<typename TestFixture::Kernel>;

		EXPECT_TRUE(approxEqual(solution[0].x(), (Vector {
			146.566682434017,
			-427.218558989821,
			-345.347969700289
		}), 1e-6));

		EXPECT_TRUE(approxEqual(solution[0].u(), (Vector {
			769.663140469139,
			-191.122524763114
		}), 1e-6));

		EXPECT_TRUE(approxEqual(solution[1].x(), (Vector {}), 1e-6));
		EXPECT_TRUE(approxEqual(solution[1].u(), (Vector {}), 1e-6));
	}

	/**
	 \brief Check that the case when some state and input bounds are infinite is correctly handled by a QP solver.
	 */
	TYPED_TEST_P(QpWorkspaceSolveTest, testInfiniteBoundsSome)
	{
		using Real = typename TestFixture::Real;

		typename TestFixture::Workspace ws {std::array<OcpSize, 2> { OcpSize(3, 2, 0), OcpSize(0, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(DynamicMatrix<typename TestFixture::Kernel>({
			{2.73449304081811,	1.88589647487061,	2.07850896302612},
			{1.88589647487061,	2.23398400982993,	2.04607087035366},
			{2.07850896302612,	2.04607087035366,	2.75909914758841}
		}));
		
		problem[0].R(DynamicMatrix<typename TestFixture::Kernel>({
			{2.58480403401811,	2.27683785836085},
			{2.27683785836085,	2.48531656389865}
		}));

		problem[0].S(DynamicMatrix<typename TestFixture::Kernel>({
			{1.94423942432824,	1.95671439063179},
			{2.31644542710892,	2.08748866537729},
			{2.46061457306194,	1.94728938270016}
		}));

		problem[0].q(DynamicVector<typename TestFixture::Kernel>({
			0.276025076998578,
			0.679702676853675,
			0.655098003973841
		}));

		problem[0].r(DynamicVector<typename TestFixture::Kernel>({
			0.162611735194631,
			0.118997681558377
		}));

		problem[0].A(StaticMatrix<typename TestFixture::Kernel, 0, 3>(0.));
		problem[0].B(StaticMatrix<typename TestFixture::Kernel, 0, 2>(0.));
		problem[0].b(StaticVector<typename TestFixture::Kernel, 0>(0.));

		// TODO: can we make lbx() etc. accept initializer lists?
		problem[0].lbx(DynamicVector<typename TestFixture::Kernel> {-10000., -inf<Real>(), -10000.});
		problem[0].ubx(DynamicVector<typename TestFixture::Kernel> {10000., inf<Real>(), 10000.});
		problem[0].lbu(DynamicVector<typename TestFixture::Kernel> {-inf<Real>(), -10000.});
		problem[0].ubu(DynamicVector<typename TestFixture::Kernel> {inf<Real>(), 10000.});

		ws.solve();
		auto solution = ws.solution();

		using Vector = DynamicVector<typename TestFixture::Kernel>;

		EXPECT_TRUE(approxEqual(solution[0].x(), (Vector {
			146.566682434017,
			-427.218558989821,
			-345.347969700289
		}), 1e-6));

		EXPECT_TRUE(approxEqual(solution[0].u(), (Vector {
			769.663140469139,
			-191.122524763114
		}), 1e-6));

		EXPECT_TRUE(approxEqual(solution[1].x(), (Vector {}), 1e-6));
		EXPECT_TRUE(approxEqual(solution[1].u(), (Vector {}), 1e-6));
	}

	TYPED_TEST_P(QpWorkspaceSolveTest, testSolve2)
	{
		using Real = typename TestFixture::Real;
		using Kernel = typename TestFixture::Kernel;
		
		typename TestFixture::Workspace workspace {OcpSize {3, 0, 0}, OcpSize {0, 0, 0}};
		
		auto& stage0 = workspace.problem()[0];
		stage0.gaussNewtonCostApproximation(
			DynamicVector<Kernel> {1., 2., 42.},
			IdentityMatrix<Kernel> {3u},
			DynamicMatrix<Kernel> {3u, 0u}
		);
		stage0.stateBounds(-inf<Real>(), inf<Real>());
		stage0.inputBounds(-inf<Real>(), inf<Real>());

		workspace.solve();
	
		EXPECT_EQ(forcePrint(workspace.solution()[0].x()), forcePrint(DynamicVector<Kernel> {-1., -2., -42.}));
	}

	REGISTER_TYPED_TEST_CASE_P(QpWorkspaceSolveTest,
		testMoveConstructor, 
		testSolve0, 
		testSolve1, 
		DISABLED_testSolve1stage1d, 
		testSolve2stage1d, 
		testSolve2stage5d, 
		testInfiniteBoundsAll, 
		testInfiniteBoundsSome,
		testSolve2
	);
}