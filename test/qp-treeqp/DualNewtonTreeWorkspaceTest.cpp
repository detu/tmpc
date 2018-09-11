#include "TreeQpWorkspaceSolveTest.hpp"
#include "TreeQpWorkspaceTest.hpp"
//#include "QpWorkspaceSolveTest.hpp"

#include <tmpc/qp/DualNewtonTreeWorkspace.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/EigenKernel.hpp>

/*
* Unit tests for common functionality of QP workspaces.
*/

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/Matrix.hpp>

#include <tmpc/test_tools.hpp>
#include <gtest/gtest.h>

#include <iostream>
#include <array>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree_Blaze_double, TreeQpWorkspaceTest, DualNewtonTreeWorkspace<BlazeKernel<double>>);
	INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree_Eigen_double, TreeQpWorkspaceTest, DualNewtonTreeWorkspace<EigenKernel<double>>);
    INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree_Blaze_double, TreeQpWorkspaceSolveTest, DualNewtonTreeWorkspace<BlazeKernel<double>>);
    INSTANTIATE_TYPED_TEST_CASE_P(DualNewtonTree_Eigen_double, TreeQpWorkspaceSolveTest, DualNewtonTreeWorkspace<EigenKernel<double>>);
}