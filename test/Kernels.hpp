// Defines LA kernels to be used in testing.
#include <tmpc/EigenKernel.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    using Kernels = ::testing::Types<
        EigenKernel<double>,
        BlazeKernel<double>
    >;
}
