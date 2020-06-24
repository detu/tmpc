#include <blaze/Math.h>

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    /// @brief Reproduces this issue:
    /// https://bitbucket.org/blaze-lib/blaze/issues/300/calling-gges-function-with-a-select
    TEST(MathTest, testGges)
    {
        static constexpr size_t N = 5;

        blaze::StaticMatrix<double, N, N> A, B, VSL, VSR;
        blaze::StaticVector<std::complex<double>, N> alpha;
        blaze::StaticVector<double, N> beta;
        
        randomize(A);
        randomize(B);

        const auto select = []( double const * alphar, double const * alphai, double const * beta ) -> int {
            return pow(*alphar, 2) + pow(*alphai, 2) > pow(*beta, 2);
        };

        gges(A, B, VSL, alpha, beta, VSR, select);
    }
}