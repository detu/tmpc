#include "QpWorkspaceTest.hpp"
#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/HpipmWorkspace.hpp>

namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(Hpipm, QpWorkspaceTest, HpipmWorkspace<double>);
	INSTANTIATE_TYPED_TEST_CASE_P(Hpipm, QpWorkspaceSolveTest, HpipmWorkspace<double>);
}