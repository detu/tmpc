#include <tmpc/system/LtiContinuousToDiscrete.hpp>
#include <tmpc/Testing.hpp>

#include <unsupported/Eigen/MatrixFunctions>

#include <iostream>


namespace tmpc :: testing
{
    TEST(LtiContinuousToDiscreteTest, testConversionDoubleIntegrator)
    {
        blaze::StaticMatrix<double, 2, 2> const Ac {
            {0., 1.},
            {0., 0.}
        };

        blaze::StaticMatrix<double, 2, 1> const Bc {
            {0.},
            {1.}
        };

        double const time_step = 0.1;

        blaze::StaticMatrix<double, 2, 2> Ad;
        blaze::StaticMatrix<double, 2, 1> Bd;

        LtiContinuousToDiscrete<double> c2d(2, 1);
        c2d(Ac, Bc, time_step, Ad, Bd);
        
        TMPC_EXPECT_EQ(Ad, (blaze::DynamicMatrix<double> {
            {1., time_step},
            {0., 1.}
        }));

        TMPC_EXPECT_EQ(Bd, (blaze::DynamicMatrix<double> {
            {pow(time_step, 2) / 2.},
            {time_step}
        }));
    }
}