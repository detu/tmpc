/// This test suite is not finished and is excluded from build.

#include "StaticRiccatiTest.hpp"

#include <tmpc/qp/StaticClassicalRiccati.hpp>
#include <tmpc/qp/StaticFactorizedRiccati.hpp>


namespace tmpc :: testing
{
    using StaticClassicalRiccati_double_2_1 = StaticClassicalRiccati<double, 2, 1>;
    using StaticFactorizedRiccati_double_2_1 = StaticFactorizedRiccati<double, 2, 1>;

    INSTANTIATE_TYPED_TEST_SUITE_P(StaticClassicalRiccati, StaticRiccatiTest, StaticClassicalRiccati_double_2_1);
    INSTANTIATE_TYPED_TEST_SUITE_P(StaticFactorizedRiccati, StaticRiccatiTest, StaticFactorizedRiccati_double_2_1);
}