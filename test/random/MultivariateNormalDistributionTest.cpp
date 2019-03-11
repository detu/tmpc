#include <tmpc/random/MultivariateNormalDistribution.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    TEST(MultivariateNormalDistributionTest, testConstructFromDimension)
    {
        size_t constexpr N = 2;
        
        MultivariateNormalDistribution<double> dist(N);
        
        EXPECT_EQ(dist.dimension(), N);
        TMPC_EXPECT_EQ(dist.covariance(), blaze::IdentityMatrix<double>(N));
        TMPC_EXPECT_EQ(dist.mean(), blaze::ZeroVector<double>(N));
    }


    TEST(MultivariateNormalDistributionTest, testConstructFromMeanAndCovariance)
    {
        size_t constexpr N = 2;
        using Vec = blaze::StaticVector<double, N, blaze::columnVector>;
        using Mat = blaze::SymmetricMatrix<blaze::StaticMatrix<double, N, N>>;

        Vec const mu {0.1, 0.2};
        Mat const sigma {
            {1., 0.5},
            {0.5, 2.}
        };
        
        MultivariateNormalDistribution<double> dist(mu, sigma);
        
        EXPECT_EQ(dist.dimension(), N);
        TMPC_EXPECT_EQ(dist.covariance(), sigma);
        TMPC_EXPECT_EQ(dist.mean(), mu);
    }


    TEST(MultivariateNormalDistributionTest, testGenerate)
    {
        size_t constexpr N = 2;
        using Vec = blaze::StaticVector<double, N, blaze::columnVector>;
        using Mat = blaze::SymmetricMatrix<blaze::StaticMatrix<double, N, N>>;

        Vec const mu {0.1, 0.2};
        Mat const sigma {
            {1., 0.5},
            {0.5, 2.}
        };
        
        MultivariateNormalDistribution<double> dist(mu, sigma);
        std::mt19937 gen;

        Vec sum_x(0.);
        Mat sum_xxT = blaze::ZeroMatrix<double>(N, N);

        size_t constexpr num_samples = 1000000;
        for (size_t i = 0; i < num_samples; ++i)
        {
            Vec const x = dist(gen);
            sum_x += x;
            sum_xxT += (x - mu) * trans(x - mu);
        }

        TMPC_EXPECT_APPROX_EQ(sum_x / num_samples, mu, 1e-3, 0.);
        TMPC_EXPECT_APPROX_EQ(sum_xxT / num_samples, sigma, 0., 1e-2);
    }
}