/// This test suite is not finished and is excluded from build.

#include "StaticRiccatiTest.hpp"

// #include <tmpc/qp/StaticClassicalRiccati.hpp>
#include <tmpc/qp/StaticFactorizedRiccati.hpp>


namespace tmpc :: testing
{
    // INSTANTIATE_TYPED_TEST_SUITE_P(ClassicalRiccati, StaticRiccatiTest, ClassicalRiccati<double>);
    INSTANTIATE_TYPED_TEST_SUITE_P(StaticFactorizedRiccati, StaticRiccatiTest, StaticFactorizedRiccati);
}