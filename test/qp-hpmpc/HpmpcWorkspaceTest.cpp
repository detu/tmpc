#include "QpWorkspaceTest.hpp"
#include "QpWorkspaceSolveTest.hpp"
#include "TreeQpWorkspaceTest.hpp"
#include "TreeQpWorkspaceSolveTest.hpp"

#include <tmpc/qp/HpmpcWorkspace.hpp>
//#include <tmpc/qp/TreeQpWorkspaceAdaptor.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/EigenKernel.hpp>

namespace tmpc :: testing
{
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Blaze_double, QpWorkspaceTest, HpmpcWorkspace<double>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Blaze_double, QpWorkspaceSolveTest, HpmpcWorkspace<double>);
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Blaze_double, TreeQpWorkspaceTest, HpmpcWorkspace<double>);
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Blaze_double, TreeQpWorkspaceSolveTest, HpmpcWorkspace<double>);
}