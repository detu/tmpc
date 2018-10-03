#include "TreeQpWorkspaceSolveTest.hpp"
#include "TreeQpWorkspaceTest.hpp"
//#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/MpipmWorkspace.hpp>

#include <tmpc/test_tools.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(Mpipm_double, TreeQpWorkspaceTest, MpipmWorkspace<double>);
	//INSTANTIATE_TYPED_TEST_CASE_P(Mpipm_double, TreeQpWorkspaceSolveTest, MpipmWorkspace<double>);

	TEST(Mpipm_double, testCtor)
	{
		size_t constexpr N = 3, NX = 2, NU = 1;
		MpipmWorkspace<double> ws(ocpGraphLinear(N + 1), ocpSizeNominalMpc(N, NX, NU));

		EXPECT_EQ(columns(get(ws.Q(), 1)), NX);
		EXPECT_EQ(rows(get(ws.Q(), 1)), NX);

		put(ws.Q(), 1, blaze::DynamicMatrix<double>(NX, NX, 42.));
	}
}