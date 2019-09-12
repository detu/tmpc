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
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Blaze_double, QpWorkspaceTest, HpmpcWorkspace<BlazeKernel<double>>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Blaze_double, QpWorkspaceSolveTest, HpmpcWorkspace<BlazeKernel<double>>);
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Blaze_double, TreeQpWorkspaceTest, HpmpcWorkspace<BlazeKernel<double>>);
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Blaze_double, TreeQpWorkspaceSolveTest, HpmpcWorkspace<BlazeKernel<double>>);

	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Eigen_double, QpWorkspaceTest, HpmpcWorkspace<EigenKernel<double>>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Eigen_double, QpWorkspaceSolveTest, HpmpcWorkspace<EigenKernel<double>>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Eigen_double, TreeQpWorkspaceTest, TreeQpWorkspaceAdaptor<HpmpcWorkspace<EigenKernel<double>>>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpmpc_Eigen_double, TreeQpWorkspaceSolveTest, TreeQpWorkspaceAdaptor<HpmpcWorkspace<EigenKernel<double>>>);
}