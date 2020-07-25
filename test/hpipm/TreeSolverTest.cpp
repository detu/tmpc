#include <test/qp/QpSolverTest.hpp>
// #include "SoftConstraintsTest.hpp"
#include "SolveUnconstrainedTest.hpp"

#include <tmpc/hpipm/TreeSolver.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_SUITE_P(TreeSolver_double, QpSolverTest, hpipm::TreeSolver<double>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(TreeSolver_double, SoftConstraintsTest, TreeSolver<double>);
	INSTANTIATE_TYPED_TEST_SUITE_P(TreeSolver_double, SolveUnconstrainedTest, hpipm::TreeSolver<double>);
}