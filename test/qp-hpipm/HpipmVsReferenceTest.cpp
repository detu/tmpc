#include <tmpc/qp/HpipmWorkspace.hpp>
#include <tmpc/qp/MpipmWorkspace.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/qp/ClassicalRiccati.hpp>
#include <tmpc/qp/FactorizedRiccati.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	// INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, TreeQpWorkspaceTest, HpipmWorkspace<BlazeKernel<double>>);
	// INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, TreeQpWorkspaceSolveTest, HpipmWorkspace<BlazeKernel<double>>);
	// //INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, SoftConstraintsTest, HpipmWorkspace<BlazeKernel<double>>);

	// INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, SolveUnconstrainedTest, HpipmWorkspace<BlazeKernel<double>>);	

	TEST(HpipmVsReferenceTest, testClassicalRiccati)
	{
		size_t const N = 5, nx = 3, nu = 2;

        OcpGraph const g = ocpGraphLinear(N + 1);
        auto const sz = ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true);
        
        HpipmWorkspace<BlazeKernel<double>> ws_hpipm(g, sz);
        MpipmWorkspace<double> ws_mpipm(g, sz);
        ClassicalRiccati<double> riccati(g, sz);

        randomizeQp(ws_mpipm);        
        copyQpProperties(ws_mpipm, ws_hpipm);

        ws_hpipm.solveUnconstrained();
        riccati(ws_mpipm, ws_mpipm);

        for (auto v : graph::vertices(g))
        {
            EXPECT_EQ(forcePrint(get(ws_mpipm.x(), v)), forcePrint(get(ws_hpipm.x(), v))) << " at node " << get(graph::vertex_index, g, v);
            EXPECT_EQ(forcePrint(get(ws_mpipm.u(), v)), forcePrint(get(ws_hpipm.u(), v))) << " at node " << get(graph::vertex_index, g, v);
        }
	}


    TEST(HpipmVsReferenceTest, testFactorizedRiccati)
	{
		size_t const N = 5, nx = 3, nu = 2;

        OcpGraph const g = ocpGraphLinear(N + 1);
        auto const sz = ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true);
        
        HpipmWorkspace<BlazeKernel<double>> ws_hpipm(g, sz);
        MpipmWorkspace<double> ws_mpipm(g, sz);
        FactorizedRiccati<double> riccati(g, sz);

        randomizeQp(ws_mpipm);        
        copyQpProperties(ws_mpipm, ws_hpipm);

        ws_hpipm.solveUnconstrained();
        riccati(ws_mpipm, ws_mpipm);

        for (auto v : graph::vertices(g))
        {
            EXPECT_EQ(forcePrint(get(ws_mpipm.x(), v)), forcePrint(get(ws_hpipm.x(), v))) << " at node " << get(graph::vertex_index, g, v);
            EXPECT_EQ(forcePrint(get(ws_mpipm.u(), v)), forcePrint(get(ws_hpipm.u(), v))) << " at node " << get(graph::vertex_index, g, v);
        }
	}
}