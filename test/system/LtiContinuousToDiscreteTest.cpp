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


    TEST(LtiContinuousToDiscreteTest, testConversion5x2)
    {
        blaze::DynamicMatrix<double> const Ac {
            {8.147236863931789e-01,     9.754040499940952e-02,     1.576130816775483e-01,     1.418863386272153e-01,     6.557406991565868e-01},
            {9.057919370756192e-01,     2.784982188670484e-01,     9.705927817606157e-01,     4.217612826262750e-01,     3.571167857418955e-02},
            {1.269868162935061e-01,     5.468815192049838e-01,     9.571669482429456e-01,     9.157355251890671e-01,     8.491293058687771e-01},
            {9.133758561390194e-01,     9.575068354342976e-01,     4.853756487228412e-01,     7.922073295595544e-01,     9.339932477575505e-01},
            {6.323592462254095e-01,     9.648885351992765e-01,     8.002804688888001e-01,     9.594924263929030e-01,     6.787351548577735e-01},
        };

        blaze::DynamicMatrix<double> const Bc {
            {7.577401305783334e-01,     7.060460880196088e-01},
            {7.431324681249162e-01,     3.183284637742068e-02},
            {3.922270195341682e-01,     2.769229849608900e-01},
            {6.554778901775566e-01,     4.617139063115394e-02},
            {1.711866878115618e-01,     9.713178123584754e-02},
        };

        double const time_step = 0.1;

        blaze::DynamicMatrix<double> Ad;
        blaze::DynamicMatrix<double> Bd;

        LtiContinuousToDiscrete<double> c2d(5, 2);
        c2d(Ac, Bc, time_step, Ad, Bd);
        
        TMPC_EXPECT_EQ(Ad, (blaze::DynamicMatrix<double> {
            {1.088741709584660e+00,     1.512689356549103e-02,     2.124624146254563e-02,     2.003599314712967e-02,     7.246034986634130e-02},
            {9.909453156694679e-02,     1.034346067893570e+00,     1.058525264843489e-01,     5.061363934412986e-02,     1.371197685629524e-02},
            {2.463671854591314e-02,     6.789324982996729e-02,     1.110170988307421e+00,     1.063538359702284e-01,     9.805188853349527e-02},
            {1.078374428042223e-01,     1.084446346078798e-01,     6.340491937361292e-02,     1.093201697172563e+00,     1.069928104736642e-01},
            {7.879902173001530e-02,     1.095350784292692e-01,     9.565259383512306e-02,     1.106801110930133e-01,     1.081830698158303e+00},
        }));

        TMPC_EXPECT_EQ(Bd, (blaze::DynamicMatrix<double> {
            {8.109899998657424e-02,     7.431769777268771e-02},
            {8.280816621049328e-02,     8.208946675009486e-03},
            {4.852287590951531e-02,     3.065557584184504e-02},
            {7.818775312987819e-02,     9.891525437807851e-03},
            {2.968117707753279e-02,     1.433470474925235e-02},
        }));
    }
}