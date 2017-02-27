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
#include <tmpc/Matrix.hpp>
#include <tmpc/core/problem_specific.hpp>

#include "qp_test_problems.hpp"
#include "gtest_tools_eigen.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <utility>

using namespace tmpc;

template <typename Solver_>
class QPSolverTest : public ::testing::Test
{
public:
	typedef Solver_ Solver;
	typedef typename Solver::Problem Problem;
	typedef typename Solver::Solution Solution;
	using Scalar = typename Solver::Scalar;

protected:

	QPSolverTest()
	:	size_(RtiQpSize(NT, NX, NU, NC, NCT))
	,	solver_(size_.begin(), size_.end())
	{		
	}

	Problem ConstructProblem()
	{
		return Problem(size_.begin(), size_.end());
	}

	Solution ConstructSolution()
	{
		return Solution(size_.begin(), size_.end());
	}

	// Define dimensions
	static unsigned constexpr NX = 2;
	static unsigned constexpr NU = 1;
	static unsigned constexpr NW = 0;
	static unsigned constexpr NY = 0;
	static unsigned constexpr NP = 0;
	static unsigned constexpr NC = 0;
	static unsigned constexpr NCT = 0;

	typedef StaticVector<Scalar, NX> StateVector;
	typedef StaticVector<Scalar, NU> InputVector;
	
	unsigned const NT = 2;
	std::vector<QpSize> size_;
	Solver solver_;
};

// Define a kernel
/*
typedef tmpc::EigenKernel<double> K;

typedef tmpc::ProblemSpecific<tmpc::EigenKernel<double>, Dimensions> PS;
*/

typedef ::testing::Types<
		tmpc::CondensingSolver<double>
,		tmpc::HPMPCSolver     <double>
,		tmpc::QpOasesSolver
	> SolverTypes;

TYPED_TEST_CASE(QPSolverTest, SolverTypes);

/// \brief Check if QPSolver move constructor works and the solver works after move constructor.
TYPED_TEST(QPSolverTest, move_constructor_test)
{
	auto const qp = tmpc_test::qp_problems::problem_0<typename TestFixture::Problem>();
	typename TestFixture::Solution solution(sizeBegin(qp), sizeEnd(qp));

	typename TestFixture::Solver solver = std::move(this->solver_);
	solver.Solve(qp, solution);

	using StateVector = typename TestFixture::StateVector;
	using InputVector = typename TestFixture::InputVector;

	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].get_x(), (StateVector {1., -1.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].get_u(), (InputVector {-1.}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[1].get_x(), (StateVector {0.5, 0.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[1].get_u(), (InputVector {-1.}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[2].get_x(), (StateVector {1., 1.}));
}

TYPED_TEST(QPSolverTest, solve_test_0)
{
	auto const qp = tmpc_test::qp_problems::problem_0<typename TestFixture::Problem>();
	typename TestFixture::Solution solution(sizeBegin(qp), sizeEnd(qp));

	this->solver_.Solve(qp, solution);

	using StateVector = typename TestFixture::StateVector;
	using InputVector = typename TestFixture::InputVector;

	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].get_x(), (StateVector {1., -1.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].get_u(), (InputVector {-1.}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[1].get_x(), (StateVector {0.5, 0.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[1].get_u(), (InputVector {-1.}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[2].get_x(), (StateVector {1., 1.}));
}

TYPED_TEST(QPSolverTest, solve_test_1)
{
	auto const qp = tmpc_test::qp_problems::problem_1<typename TestFixture::Problem>();
	typename TestFixture::Solution solution(sizeBegin(qp), sizeEnd(qp));

	this->solver_.Solve(qp, solution);

	using StateVector = typename TestFixture::StateVector;
	using InputVector = typename TestFixture::InputVector;

	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].get_x(), (StateVector {1., 0.}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[0].get_u(), (InputVector {-0.690877362606266}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[1].get_x(), (StateVector {0.654561318696867, -0.690877362606266}));
	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[1].get_u(), (InputVector {0.215679569867116}));

	EXPECT_PRED2(MatrixApproxEquality(1e-6), solution[2].get_x(), (StateVector {0.0715237410241597, -0.475197792739149}));
}
