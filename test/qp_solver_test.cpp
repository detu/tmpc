/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include <tmpc/qp/CondensingSolver.hpp>
#include <tmpc/qp/HPMPCSolver.hpp>
#include <tmpc/qp/QpOasesSolver.hpp>
#include <tmpc/core/RealtimeIteration.hpp>	// for RtiQpSize(0
#include <tmpc/qp/Printing.hpp>
#include <tmpc/kernel/eigen.hpp>
#include <tmpc/core/problem_specific.hpp>

#include "qp_test_problems.hpp"
#include "gtest_tools_eigen.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <utility>

template <typename Solver_>
class QPSolverTest : public ::testing::Test
{
public:
	typedef Solver_ Solver;
	typedef typename Solver::Problem Problem;
	typedef typename Solver::Solution Solution;

protected:

	QPSolverTest();
	Problem ConstructProblem();
	Solution ConstructSolution();

	unsigned const NT = 2;
	Solver solver_;
};

template <typename Solver_>
QPSolverTest<Solver_>::QPSolverTest()
:	solver_(NT)
{
}

template <typename Solver_>
typename QPSolverTest<Solver_>::Problem QPSolverTest<Solver_>::ConstructProblem()
{
	return Problem(NT);
}

template <typename Solver_>
typename QPSolverTest<Solver_>::Solution QPSolverTest<Solver_>::ConstructSolution()
{
	return Solution(NT);
}

// Define dimensions
struct Dimensions
{
	static unsigned constexpr NX = 2;
	static unsigned constexpr NU = 1;
	static unsigned constexpr NW = 0;
	static unsigned constexpr NY = 0;
	static unsigned constexpr NP = 0;
	static unsigned constexpr NC = 0;
	static unsigned constexpr NCT = 0;
};

template <>
QPSolverTest<tmpc::QpOasesSolver>::QPSolverTest()
:	solver_(tmpc::RtiQpSize(NT, Dimensions::NX, Dimensions::NU, Dimensions::NC, Dimensions::NCT))
{
}

template <>
typename QPSolverTest<tmpc::QpOasesSolver>::Problem QPSolverTest<tmpc::QpOasesSolver>::ConstructProblem()
{
	return Problem(solver_.size());
}

template <>
typename QPSolverTest<tmpc::QpOasesSolver>::Solution QPSolverTest<tmpc::QpOasesSolver>::ConstructSolution()
{
	return Solution(solver_.size());
}

// Define a kernel
typedef tmpc::EigenKernel<double> K;

typedef tmpc::ProblemSpecific<tmpc::EigenKernel<double>, Dimensions> PS;

typedef ::testing::Types<
		tmpc::CondensingSolver<K, Dimensions>
,		tmpc::HPMPCSolver     <K, Dimensions>
,		tmpc::QpOasesSolver
	> SolverTypes;

TYPED_TEST_CASE(QPSolverTest, SolverTypes);

/// \brief Check if QPSolver move constructor works and the solver works after move constructor.
TYPED_TEST(QPSolverTest, move_constructor_test)
{
	auto qp = this->ConstructProblem();
	tmpc_test::qp_problems::problem_0(qp);

	auto solution = this->ConstructSolution();

	typename TestFixture::Solver solver = std::move(this->solver_);
	solver.Solve(qp, solution);

	PS::StateVector x_expected;
	PS::InputVector u_expected;

	x_expected << 1., -1.;	u_expected << -1;
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_x(0), x_expected);
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_u(0), u_expected);

	x_expected << 0.5, 0.;	u_expected << -1;
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_x(1), x_expected);
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_u(1), u_expected);

	x_expected << 1., 1;
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_x(2), x_expected);
}

TYPED_TEST(QPSolverTest, solve_test_0)
{
	auto qp = this->ConstructProblem();
	tmpc_test::qp_problems::problem_0(qp);

	auto solution = this->ConstructSolution();
	this->solver_.Solve(qp, solution);

	PS::StateVector x_expected;
	PS::InputVector u_expected;

	x_expected << 1., -1.;	u_expected << -1;
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_x(0), x_expected);
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_u(0), u_expected);

	x_expected << 0.5, 0.;	u_expected << -1;
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_x(1), x_expected);
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_u(1), u_expected);

	x_expected << 1., 1;
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_x(2), x_expected);
}

TYPED_TEST(QPSolverTest, solve_test_1)
{
	auto qp = this->ConstructProblem();
	tmpc_test::qp_problems::problem_1(qp);

	auto solution = this->ConstructSolution();
	this->solver_.Solve(qp, solution);

	PS::StateVector x_expected;
	PS::InputVector u_expected;

	x_expected << 1., 0.;	u_expected << -0.690877362606266;
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_x(0), x_expected);
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_u(0), u_expected);

	x_expected << 0.654561318696867, -0.690877362606266;	u_expected << 0.215679569867116;
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_x(1), x_expected);
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_u(1), u_expected);

	x_expected << 0.0715237410241597, -0.475197792739149;
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution.get_x(2), x_expected);
}
