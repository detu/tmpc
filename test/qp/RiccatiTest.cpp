#include "RiccatiTest.hpp"

#include <tmpc/qp/ClassicalRiccati.hpp>


namespace tmpc :: testing
{
    INSTANTIATE_TYPED_TEST_CASE_P(ClassicalRiccati, RiccatiTest, ClassicalRiccati<double>);
}