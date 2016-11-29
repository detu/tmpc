/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include <tmpc/qp/CondensingSolver.hpp>
#include <tmpc/qp/HPMPCSolver.hpp>
#include <tmpc/qp/Printing.hpp>
#include <tmpc/kernel/eigen.hpp>
#include <tmpc/core/problem_specific.hpp>

#include "qp_test_problems.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <utility>

template <typename Solver_>
class QPSolverTest : public ::testing::Test
{
public:
	typedef Solver_ Solver;

protected:
	unsigned const NT = 2;

	QPSolverTest()
	:	solver_(NT)
	{
	}

	Solver solver_;
};

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

// Define a kernel
typedef tmpc::EigenKernel<double> K;

typedef tmpc::ProblemSpecific<tmpc::EigenKernel<double>, Dimensions> PS;

typedef ::testing::Types<
		tmpc::CondensingSolver<K, Dimensions>
,		tmpc::HPMPCSolver     <K, Dimensions>
	> SolverTypes;

TYPED_TEST_CASE(QPSolverTest, SolverTypes);

/// \brief Check if QPSolver move constructor works and the solver works after move constructor.
TYPED_TEST(QPSolverTest, move_constructor_test)
{
	typename TestFixture::Solver::Problem qp(this->solver_.nT());
	tmpc_test::qp_problems::problem_0(qp);

	typename TestFixture::Solver::Solution solution(this->solver_.nT());

	typename TestFixture::Solver solver = std::move(this->solver_);
	solver.Solve(qp, solution);

	PS::StateVector x_expected;
	PS::InputVector u_expected;

	x_expected << 1., -1.;	u_expected << -1;
	EXPECT_TRUE(solution.get_x(0).isApprox(x_expected));
	EXPECT_TRUE(solution.get_u(0).isApprox(u_expected));

	x_expected << 0.5, 0.;	u_expected << -1;
	EXPECT_TRUE(solution.get_x(1).isApprox(x_expected));
	EXPECT_TRUE(solution.get_u(1).isApprox(u_expected));

	x_expected << 1., 1;
	EXPECT_TRUE(solution.get_x(2).isApprox(x_expected));
}

TYPED_TEST(QPSolverTest, solve_test_0)
{
	typename TestFixture::Solver::Problem qp(this->solver_.nT());
	tmpc_test::qp_problems::problem_0(qp);

	typename TestFixture::Solver::Solution solution(this->solver_.nT());
	this->solver_.Solve(qp, solution);

	PS::StateVector x_expected;
	PS::InputVector u_expected;

	x_expected << 1., -1.;	u_expected << -1;
	EXPECT_TRUE(solution.get_x(0).isApprox(x_expected));
	EXPECT_TRUE(solution.get_u(0).isApprox(u_expected));

	x_expected << 0.5, 0.;	u_expected << -1;
	EXPECT_TRUE(solution.get_x(1).isApprox(x_expected));
	EXPECT_TRUE(solution.get_u(1).isApprox(u_expected));

	x_expected << 1., 1;
	EXPECT_TRUE(solution.get_x(2).isApprox(x_expected));
}

TYPED_TEST(QPSolverTest, solve_test_1)
{
	typename TestFixture::Solver::Problem qp(this->solver_.nT());
	tmpc_test::qp_problems::problem_1(qp);

	typename TestFixture::Solver::Solution solution(this->solver_.nT());
	this->solver_.Solve(qp, solution);

	PS::StateVector x_expected;
	PS::InputVector u_expected;

	x_expected << 1., 0.;	u_expected << -0.690877362606266;
	EXPECT_TRUE(solution.get_x(0).isApprox(x_expected, 1e-6));
	EXPECT_TRUE(solution.get_u(0).isApprox(u_expected, 1e-6));

	x_expected << 0.654561318696867, -0.690877362606266;	u_expected << 0.215679569867116;
	EXPECT_TRUE(solution.get_x(1).isApprox(x_expected, 1e-6));
	EXPECT_TRUE(solution.get_u(1).isApprox(u_expected, 1e-6));

	x_expected << 0.0715237410241597, -0.475197792739149;
	EXPECT_TRUE(solution.get_x(2).isApprox(x_expected, 1e-6));
}
