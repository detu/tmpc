#include <test/Kernels.hpp>

namespace tmpc :: testing
{
    template <typename Kernel_>
    class KernelTest 
    :   public Test
    {
    protected:
        using Kernel = Kernel_;
    };

    TYPED_TEST_CASE(KernelTest, Kernels);
}
