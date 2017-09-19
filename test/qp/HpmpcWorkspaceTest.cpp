#include "QpWorkspaceTest.hpp"
#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/HpmpcWorkspace.hpp>
#include <tmpc/BlazeKernel.hpp>

namespace tmpc :: testing
{
	using Kernel = BlazeKernel<double>;
	INSTANTIATE_TYPED_TEST_CASE_P(Hpmpc, QpWorkspaceTest, HpmpcWorkspace<Kernel>);
	INSTANTIATE_TYPED_TEST_CASE_P(Hpmpc, QpWorkspaceSolveTest, HpmpcWorkspace<Kernel>);
}