#include "DynamicRiccatiTest.hpp"

#include <tmpc/qp/DynamicClassicalRiccati.hpp>
#include <tmpc/qp/DynamicFactorizedRiccati.hpp>


namespace tmpc :: testing
{
    INSTANTIATE_TYPED_TEST_SUITE_P(DynamicClassicalRiccati, DynamicRiccatiTest, DynamicClassicalRiccati<double>);
    INSTANTIATE_TYPED_TEST_SUITE_P(DynamicFactorizedRiccati, DynamicRiccatiTest, DynamicFactorizedRiccati<double>);
}