#include "CondensingTest.hpp"

#include <tmpc/qp/CondensingN2.hpp>

#include <test/Kernels.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(CondensingN2_Eigen_double, CondensingTest, CondensingN2<EigenKernel<double>>);
	INSTANTIATE_TYPED_TEST_CASE_P(CondensingN2_Blaze_double, CondensingTest, CondensingN2<BlazeKernel<double>>);
}