#include "QpSolverTest.hpp"

#include <tmpc/qp/HpmpcSolver.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_double, QpSolverTest, HpmpcSolver<double>);
}
