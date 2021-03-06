#include "QpSolverTest.hpp"
//#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/DualNewtonTreeWorkspace.hpp>
#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_SUITE_P(DualNewtonTree, QpSolverTest, DualNewtonTreeWorkspace);
}