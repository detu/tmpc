#include "QpWorkspaceTest.hpp"
#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/HpmpcWorkspace.hpp>

namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(Hpmpc, QpWorkspaceTest, HpmpcWorkspace<double>);
	INSTANTIATE_TYPED_TEST_CASE_P(Hpmpc, QpWorkspaceSolveTest, HpmpcWorkspace<double>);
}