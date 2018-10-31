#include "TreeQpWorkspaceSolveTest.hpp"
#include "TreeQpWorkspaceTest.hpp"
//#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/DualNewtonTreeWorkspace.hpp>
#include <tmpc/ocp/OcpGraph.hpp>

#include <tmpc/test_tools.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree, TreeQpWorkspaceTest, DualNewtonTreeWorkspace);
    INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree, TreeQpWorkspaceSolveTest, DualNewtonTreeWorkspace);
}