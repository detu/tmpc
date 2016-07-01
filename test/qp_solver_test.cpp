/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include "../include/qp/CondensingSolver.hpp"
#include "../include/qp/HPMPCSolver.hpp"
#include "../include/qp/Printing.hpp"

#include "qp_test_problems.hpp"

#include <gtest/gtest.h>

#include <iostream>

template <typename Solver_>
class QPSolverTest : public ::testing::Test
{
public:
	typedef Solver_ Solver;
	//typedef typename Problem::StateInputVector StateInputVector;
	//typedef typename Problem::StageHessianMatrix StageHessianMatrix;
	//typedef typename Problem::InterStageMatrix InterStageMatrix;

protected:
	/*
	unsigned const NX = n_x<QP>();
	unsigned const NU = n_u<QP>();
	unsigned const NC = n_d<QP>();
	unsigned const NCT = n_d_end<QP>();
	*/

	/*
	typedef Eigen::Matrix<double, QP::NX, 1> StateVector;
	typedef Eigen::Matrix<double, QP::NU, 1> InputVector;
	typedef Eigen::Matrix<double, QP::NX, QP::NX> StateStateMatrix;
	typedef Eigen::Matrix<double, QP::NX, QP::NU> StateInputMatrix;
	typedef Eigen::Matrix<double, QP::NU, QP::NU> InputInputMatrix;
	*/

	unsigned const NT = 2;

	QPSolverTest()
	:	solver_(NT)
	{
	}

	Solver solver_;
};

typedef ::testing::Types<
		tmpc::CondensingSolver<2, 1, 0, 0>
,		tmpc::HPMPCSolver     <2, 1, 0, 0>
	> SolverTypes;

TYPED_TEST_CASE(QPSolverTest, SolverTypes);

TYPED_TEST(QPSolverTest, solve_test_0)
{
	typename TestFixture::Solver::Problem qp(this->solver_.nT());
	tmpc_test::qp_problems::problem_0(qp);

	typename TestFixture::Solver::Solution solution(this->solver_.nT());
	this->solver_.Solve(qp, solution);

	Eigen::Matrix<double, TestFixture::Solver::NX, 1> x_expected;
	Eigen::Matrix<double, TestFixture::Solver::NU, 1> u_expected;

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

	Eigen::Matrix<double, TestFixture::Solver::NX, 1> x_expected;
	Eigen::Matrix<double, TestFixture::Solver::NU, 1> u_expected;

	x_expected << 1., 0.;	u_expected << -0.690877362606266;
	EXPECT_TRUE(solution.get_x(0).isApprox(x_expected, 1e-6));
	EXPECT_TRUE(solution.get_u(0).isApprox(u_expected, 1e-6));

	x_expected << 0.654561318696867, -0.690877362606266;	u_expected << 0.215679569867116;
	EXPECT_TRUE(solution.get_x(1).isApprox(x_expected, 1e-6));
	EXPECT_TRUE(solution.get_u(1).isApprox(u_expected, 1e-6));

	x_expected << 0.0715237410241597, -0.475197792739149;
	EXPECT_TRUE(solution.get_x(2).isApprox(x_expected, 1e-6));
}
