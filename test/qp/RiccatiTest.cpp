#include "RiccatiTest.hpp"

#include <tmpc/qp/ClassicalRiccati.hpp>
#include <tmpc/qp/FactorizedRiccati.hpp>


namespace tmpc :: testing
{
    INSTANTIATE_TYPED_TEST_SUITE_P(ClassicalRiccati, RiccatiTest, ClassicalRiccati<double>);
    INSTANTIATE_TYPED_TEST_SUITE_P(FactorizedRiccati, RiccatiTest, FactorizedRiccati<double>);
}