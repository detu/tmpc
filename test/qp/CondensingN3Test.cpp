#include "CondensingTest.hpp"

#include <tmpc/qp/CondensingN3.hpp>

#include <test/Kernels.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_SUITE_P(CondensingN3_Eigen_double, CondensingTest, CondensingN3<EigenKernel<double>>);
	INSTANTIATE_TYPED_TEST_SUITE_P(CondensingN3_Blaze_double, CondensingTest, CondensingN3<BlazeKernel<double>>);
}