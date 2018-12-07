#include <tmpc/qp/HpipmWorkspace.hpp>
#include <tmpc/qp/MpipmWorkspace.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/qp/ClassicalRiccati.hpp>
#include <tmpc/qp/OcpQp.hpp>

#include <tmpc/BlazeKernel.hpp>

#include <tmpc/test_tools.hpp>


namespace tmpc :: testing
{
	// INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, TreeQpWorkspaceTest, HpipmWorkspace<BlazeKernel<double>>);
	// INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, TreeQpWorkspaceSolveTest, HpipmWorkspace<BlazeKernel<double>>);
	// //INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, SoftConstraintsTest, HpipmWorkspace<BlazeKernel<double>>);

	// INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, SolveUnconstrainedTest, HpipmWorkspace<BlazeKernel<double>>);	

	TEST(HpipmVsReferenceTest, testRiccati)
	{
		size_t const N = 1, nx = 2, nu = 1;

        OcpGraph const g = ocpGraphLinear(N + 1);
        auto const sz = ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true);
        
        HpipmWorkspace<BlazeKernel<double>> ws_hpipm(g, sz);
        MpipmWorkspace<double> ws_mpipm(g, sz);
        ClassicalRiccati<double> riccati(g, sz);

        randomizeQp(ws_hpipm);

        // Disable openblas multithreading
        openblas_set_num_threads(1);

        ws_hpipm.solveUnconstrained();
        riccati(ws_hpipm, ws_mpipm);

        for (auto v : graph::vertices(g))
        {
            EXPECT_EQ(forcePrint(get(ws_mpipm.x(), v)), forcePrint(get(ws_hpipm.x(), v)));
        }
	}
}