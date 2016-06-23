#include <qp/CondensingSolver.hpp>
#include <qp/Condensing.hpp>
#include <qp/qpOASESProgram.hpp>

#include "qp_test_problems.hpp"

#include <gtest/gtest.h>

#include <iostream>

unsigned const NX = 2;
unsigned const NU = 1;
unsigned const NC = 0;
unsigned const NCT = 0;
unsigned const NT = 2;

typedef tmpc::CondensingSolver<NX, NU, NC, NCT> Solver;
typedef Solver::Problem Problem;
typedef Solver::Solution Solution;

std::ostream& operator<<(std::ostream& os, Solution const& point)
{
	//typedef typename camels::CondensingSolver<NX_, NU_, NC_, NCT_>::size_type size_type;
	typedef unsigned size_type;
	for (size_type i = 0; i < point.nT(); ++i)
		os << point.w(i) << std::endl;

	return os << point.wend() << std::endl;
}

TEST(CondensingSolver_test, condensing_test)
{
	Solver::Problem qp(NT);
	qp.zMin(0)  .setConstant(-1);	qp.zMax(0)  .setConstant(1);
	qp.zMin(1)  .setConstant(-1);	qp.zMax(1)  .setConstant(1);
	qp.zendMin().setConstant(-1);	qp.zendMax().setConstant(1);

	// Stage 0
	Eigen::MatrixXd H0(qp.nZ(), qp.nZ());
	H0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H0 = H0.transpose() * H0;	// Make positive definite.

	const Eigen::MatrixXd Q0 = H0.topLeftCorner(qp.nX(), qp.nX());
	const Eigen::MatrixXd R0 = H0.bottomRightCorner(qp.nU(), qp.nU());
	const Eigen::MatrixXd S0 = H0.topRightCorner(qp.nX(), qp.nU());
	const Eigen::MatrixXd S0T = H0.bottomLeftCorner(qp.nU(), qp.nX());

	Eigen::MatrixXd A0(qp.nX(), qp.nX());
	A0 << 1, 1, 0, 1;

	Eigen::MatrixXd B0(qp.nX(), qp.nU());
	B0 << 0.5, 1.0;

	Eigen::VectorXd a0(qp.nX());
	a0 << 1, 2;

	// Stage 1
	Eigen::MatrixXd H1(qp.nZ(), qp.nZ());
	H1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	H1 = H1.transpose() * H1;	// Make positive definite.

	const Eigen::MatrixXd Q1 = H1.topLeftCorner(qp.nX(), qp.nX());
	const Eigen::MatrixXd R1 = H1.bottomRightCorner(qp.nU(), qp.nU());
	const Eigen::MatrixXd S1 = H1.topRightCorner(qp.nX(), qp.nU());
	const Eigen::MatrixXd S1T = H1.bottomLeftCorner(qp.nU(), qp.nX());

	Eigen::MatrixXd A1(qp.nX(), qp.nX());
	A1 << 1, 1, 0, 1;

	Eigen::MatrixXd B1(qp.nX(), qp.nU());
	B1 << 0.5, 1.0;

	Eigen::VectorXd a1(qp.nX());
	a1 << 1, 2;

	// Stage 2
	Eigen::MatrixXd H2(qp.nX(), qp.nX());
	H2 << 1, 2, 3, 4;
	H2 = H2.transpose() * H2;	// Make positive definite.

	const Eigen::MatrixXd Q2 = H2.topLeftCorner(qp.nX(), qp.nX());

	// Setup QP
	qp.H(0) = H0;
	qp.H(1) = H1;
	qp.Hend() = H2;

	qp.C(0) << A0, B0;
	qp.c(0) << a0;
	qp.C(1) << A1, B1;
	qp.c(1) << a1;

	// Condense
	camels::qpOASESProgram condensed(nIndep(qp), nDep(qp) + nConstr(qp));
	camels::Condense(qp, condensed);

	const auto Hc = condensed.H();
	Eigen::MatrixXd Hc_expected(nIndep(qp), nIndep(qp));
	
	Hc_expected <<
		A0.transpose() * Q1 * A0 + A0.transpose() * A1.transpose() * Q2 * A1 * A0 + Q0,				A0.transpose() * Q1 * B0 + A0.transpose() * A1.transpose() * Q2 * A1 * B0 + S0,	A0.transpose() * S1 + A0.transpose() * A1.transpose() * Q2 * B1,
		B0.transpose() * Q1 * A0 + B0.transpose() * A1.transpose() * Q2 * A1 * A0 + S0.transpose(), B0.transpose() * Q1 * B0 + B0.transpose() * A1.transpose() * Q2 * A1 * B0 + R0, B0.transpose() * S1 + B0.transpose() * A1.transpose() * Q2 * B1,
		S1.transpose() * A0 + B1.transpose() * Q2 * A1 * A0,										S1.transpose() * B0 + B1.transpose() * Q2 * A1 * B0,							B1.transpose() * Q2 * B1 + R1;

	EXPECT_TRUE(Hc_expected == Hc);
	std::cout << "****** Condensed problem *******" << std::endl;
	Print_MATLAB(std::cout, condensed, "qp");
	std::cout << "********* Hc_expected **********" << std::endl;
	std::cout << Hc_expected << std::endl;
	//qp.PrintQP_C(std::cout);

	const auto gc = condensed.g();
	Eigen::VectorXd gc_expected(nIndep(qp));
	gc_expected <<
		A0.transpose() * Q1 * a0					+ A0.transpose() * A1.transpose() * Q2 * a1			+ A0.transpose() * A1.transpose() * Q2 * A1 * a0,
		B0.transpose() * A1.transpose() * Q2 * a1	+ B0.transpose() * A1.transpose() * Q2 * A1 * a0	+ B0.transpose() * Q1 * a0,
		B1.transpose() * Q2 * A1 * a0				+ B1.transpose() * Q2 * a1							+ S1.transpose() * a0;

	EXPECT_TRUE(gc_expected == gc);

	/*
	std::cout << "--- A ---" << std::endl << condensed.A() << std::endl;
	std::cout << "-- lbA --" << std::endl << condensed.lbA() << std::endl;
	std::cout << "-- ubA --" << std::endl << condensed.ubA() << std::endl;
	std::cout << "-- lb ---" << std::endl << condensed.lb() << std::endl;
	std::cout << "-- ub ---" << std::endl << condensed.ub() << std::endl;
	*/
}

TEST(CondensingSolver_test, Solve_test_0)
{
	Problem qp(NT);
	tmpc_test::qp_problems::problem_0(qp);

	Solver solver(qp.nT());
	Solution solution(solver.nT());

	try
	{
		solver.Solve(qp, solution);
	}
	catch(Solver::SolveException const& x)
	{
		std::cerr << "+++++++ Condensed QP that failed: ++++++++" << std::endl;
		Print_MATLAB(std::cerr, x.getCondensedQP(), "qp");
		throw;
	}

	Solution::StateInputVector z0_expected;
	z0_expected << 1., -1., -1;
	EXPECT_TRUE(solution.w(0).isApprox(z0_expected));

	Solution::StateInputVector z1_expected;
	z1_expected << 0.5, 0., -1;
	EXPECT_TRUE(solution.w(1).isApprox(z1_expected));

	Solution::StateVector z2_expected;
	z2_expected << 1., 1;
	EXPECT_TRUE(solution.wend().isApprox(z2_expected));

	std::cout << "-- sol (condensed ) --" << std::endl << solver.getCondensedSolution() << std::endl;
	std::cout << "-- sol (multistage) --" << std::endl << solution << std::endl;
}

TEST(CondensingSolver_test, Solve_test_1)
{
	Problem qp(NT);
	tmpc_test::qp_problems::problem_1(qp);

	Solver solver(qp.nT());
	Solution solution(solver.nT());

	try
	{
		solver.Solve(qp, solution);
	}
	catch(Solver::SolveException const& x)
	{
		std::cerr << "+++++++ Condensed QP that failed: ++++++++" << std::endl;
		Print_MATLAB(std::cerr, x.getCondensedQP(), "qp");
		throw;
	}

	Solution::StateInputVector z0_expected;
	z0_expected << 1., 0., -0.690877362606266;
	EXPECT_TRUE(solution.w(0).isApprox(z0_expected));

	Solution::StateInputVector z1_expected;
	z1_expected << 0.654561318696867, -0.690877362606266, 0.215679569867116;
	EXPECT_TRUE(solution.w(1).isApprox(z1_expected));

	Solution::StateVector z2_expected;
	z2_expected << 0.0715237410241597, -0.475197792739149;
	EXPECT_TRUE(solution.wend().isApprox(z2_expected));

	std::cout << "-- sol (condensed ) --" << std::endl << solver.getCondensedSolution() << std::endl;
	std::cout << "-- sol (multistage) --" << std::endl << solution << std::endl;
}
