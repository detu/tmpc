#include "QpWorkspaceTest.hpp"
#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/HpipmWorkspace.hpp>
#include <tmpc/EigenKernel.hpp>
#include <tmpc/BlazeKernel.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Eigen_double, QpWorkspaceTest, HpipmWorkspace<EigenKernel<double>>);
	INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Eigen_double, QpWorkspaceSolveTest, HpipmWorkspace<EigenKernel<double>>);

	INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, QpWorkspaceTest, HpipmWorkspace<BlazeKernel<double>>);
	INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, QpWorkspaceSolveTest, HpipmWorkspace<BlazeKernel<double>>);
}