#include "TreeQpWorkspaceTest.hpp"
#include "TreeQpWorkspaceSolveTest.hpp"
#include "SoftConstraintsTest.hpp"
#include "SolveUnconstrainedTest.hpp"

#include <tmpc/qp/HpipmWorkspace.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpipm_double, TreeQpWorkspaceTest, HpipmWorkspace<double>);
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpipm_double, TreeQpWorkspaceSolveTest, HpipmWorkspace<double>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(Hpipm_Blaze_double, SoftConstraintsTest, HpipmWorkspace<double>);
	INSTANTIATE_TYPED_TEST_SUITE_P(Hpipm_double, SolveUnconstrainedTest, HpipmWorkspace<double>);	
}