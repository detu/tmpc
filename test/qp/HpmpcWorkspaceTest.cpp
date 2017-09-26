#include "QpWorkspaceTest.hpp"
#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/HpmpcWorkspace.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/EigenKernel.hpp>

namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(Hpmpc_Blaze_double, QpWorkspaceTest, HpmpcWorkspace<BlazeKernel<double>>);
	INSTANTIATE_TYPED_TEST_CASE_P(Hpmpc_Blaze_double, QpWorkspaceSolveTest, HpmpcWorkspace<BlazeKernel<double>>);

	INSTANTIATE_TYPED_TEST_CASE_P(Hpmpc_Eigen_double, QpWorkspaceTest, HpmpcWorkspace<EigenKernel<double>>);
	INSTANTIATE_TYPED_TEST_CASE_P(Hpmpc_Eigen_double, QpWorkspaceSolveTest, HpmpcWorkspace<EigenKernel<double>>);
}