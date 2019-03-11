#include <tmpc/numeric/CentralDifferenceDerivative.hpp>
#include <tmpc/Testing.hpp>

#include <cmath>


namespace tmpc :: testing
{
    TEST(CentralDifferenceDerivativeTest, testDiff)
    {
        size_t constexpr NX = 2;
        size_t constexpr NY = 3;

        using VecX = blaze::StaticVector<double, NX, blaze::columnVector>;
        using VecY = blaze::StaticVector<double, NY, blaze::columnVector>;

        // Test function
        auto const f = [] (VecX const& x)
        {
            return VecY {
                pow(x[0], 2) + 2. * x[0] * x[1] + 3. * x[1],
                cos(3. * x[0]) + sin(x[1]) + exp(0.1 * x[0] + 0.2 * x[1]),
                x[0] / x[1]
            };
        };

        // Analytic Jacobian of f
        auto const Jf = [] (VecX const& x)
        {
            return blaze::StaticMatrix<double, NY, NX> {
                {2. * x[0] + 2. * x[1], 2. * x[0] + 3.},
                {-3. * sin(3. * x[0]) + 0.1 * exp(0.1 * x[0] + 0.2 * x[1]), cos(x[1]) + 0.2 * exp(0.1 * x[0] + 0.2 * x[1])},
                {1. / x[1], -x[0] / pow(x[1], 2)}
            };
        };

        // Numerical differentiator
        CentralDifferenceDerivative<double> diff(NX, NY);

        // Test point
        VecX const x {1., 2.};

        // Difference delta
        VecX const delta {1e-6, 1e-6};

        // You will be surprised, by they are really equal within machine precision.
        EXPECT_EQ(forcePrint(diff(f, x, delta)), forcePrint(Jf(x)));
    }
}