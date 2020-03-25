#include "QpWorkspaceTest.hpp"
#include "QpWorkspaceSolveTest.hpp"
#include "TreeQpWorkspaceTest.hpp"
#include "TreeQpWorkspaceSolveTest.hpp"

#include <tmpc/qp/HpmpcWorkspace.hpp>
//#include <tmpc/qp/TreeQpWorkspaceAdaptor.hpp>


namespace tmpc :: testing
{
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_double, QpWorkspaceTest, HpmpcWorkspace<double>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_double, QpWorkspaceSolveTest, HpmpcWorkspace<double>);
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_double, TreeQpWorkspaceTest, HpmpcWorkspace<double>);
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_double, TreeQpWorkspaceSolveTest, HpmpcWorkspace<double>);
}