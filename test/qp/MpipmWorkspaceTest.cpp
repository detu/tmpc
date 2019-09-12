#include "TreeQpWorkspaceSolveTest.hpp"
#include "TreeQpWorkspaceTest.hpp"

#include <tmpc/qp/MpipmWorkspace.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_SUITE_P(Mpipm_double, TreeQpWorkspaceTest, MpipmWorkspace<double>);
	//INSTANTIATE_TYPED_TEST_SUITE_P(Mpipm_double, TreeQpWorkspaceSolveTest, MpipmWorkspace<double>);
}