#include <test/qp/QpSolverTest.hpp>
// #include "SoftConstraintsTest.hpp"
#include "SolveUnconstrainedTest.hpp"

#include <tmpc/hpipm/NominalSolver.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_SUITE_P(NominalSolver_double, QpSolverTest, hpipm::NominalSolver<double>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(NominalSolver_double, SoftConstraintsTest, NominalSolver<double>);
	INSTANTIATE_TYPED_TEST_SUITE_P(NominalSolver_double, SolveUnconstrainedTest, hpipm::NominalSolver<double>);

	// INSTANTIATE_TYPED_TEST_SUITE_P(NominalSolver_float, QpSolverTest, NominalSolver<float>);
	// INSTANTIATE_TYPED_TEST_SUITE_P(NominalSolver_float, TreeQpSolverTest, NominalSolver<float>);
	// //INSTANTIATE_TYPED_TEST_SUITE_P(NominalSolver_float, SoftConstraintsTest, NominalSolver<float>);
	// INSTANTIATE_TYPED_TEST_SUITE_P(NominalSolver_float, SolveUnconstrainedTest, NominalSolver<float>);
}