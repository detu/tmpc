#include "TreeQpWorkspaceSolveTest.hpp"
#include "TreeQpWorkspaceTest.hpp"
//#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/DualNewtonTreeWorkspace.hpp>
#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_SUITE_P(DualNewtonTree, TreeQpWorkspaceTest, DualNewtonTreeWorkspace);
    INSTANTIATE_TYPED_TEST_SUITE_P(DualNewtonTree, TreeQpWorkspaceSolveTest, DualNewtonTreeWorkspace);
}