#include "TreeQpWorkspaceTest.hpp"
#include "TreeQpWorkspaceSolveTest.hpp"
#include "SoftConstraintsTest.hpp"
#include "../qp/SolveUnconstrainedTest.hpp"

#include <tmpc/qp/HpipmWorkspace.hpp>
//#include <tmpc/EigenKernel.hpp>
#include <tmpc/BlazeKernel.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpipm_Blaze_double, TreeQpWorkspaceTest, HpipmWorkspace<BlazeKernel<double>>);
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpipm_Blaze_double, TreeQpWorkspaceSolveTest, HpipmWorkspace<BlazeKernel<double>>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpipm_Blaze_double, SoftConstraintsTest, HpipmWorkspace<BlazeKernel<double>>);

	INSTANTIATE_TYPED_TEST_SUITE_P(Hpipm_Blaze_double, SolveUnconstrainedTest, HpipmWorkspace<BlazeKernel<double>>);	
}