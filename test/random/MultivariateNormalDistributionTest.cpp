#include <tmpc/random/MultivariateNormalDistribution.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    TEST(MultivariateNormalDistributionTest, test)
    {
        size_t constexpr N = 2;
        
        MultivariateNormalDistribution<double> dist(N);
        std::mt19937 gen;

        dist(gen);
    }
}