#include <tmpc/estimation/KalmanFilter.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    using Real = double;
    using Mat = blaze::DynamicMatrix<Real>;
    using Vec = blaze::DynamicVector<Real, blaze::columnVector>;
    
    using blaze::columnVector;
                
    template <size_t N, bool TF = blaze::columnVector>
    using StaticVector = blaze::StaticVector<Real, N, TF>;

    template <size_t M, size_t N, bool SO = blaze::columnMajor>
    using StaticMatrix = blaze::StaticMatrix<Real, M, N, SO>;


    class KalmanFilterInitTest
    :   public Test
    {
    protected:
        static size_t constexpr NX = 4;
        static size_t constexpr NY = 2;

        KalmanFilter<Real> kalman_ {NX, NY};
    };


    class KalmanFilterTest
    :   public Test
    {
    protected:
        KalmanFilterTest()
        :   kalman_ {NX, NY}
        ,   x0_ {0.8244,    0.9827,    0.7302,    0.3439}
        ,   u0_ {0.5841,    0.1078,    0.9063}
        ,   A_ {
                {0.1734,    0.0605,    0.6569,    0.0155},
                {0.3909,    0.3993,    0.6280,    0.9841},
                {0.8314,    0.5269,    0.2920,    0.1672},
                {0.8034,    0.4168,    0.4317,    0.1062}
            }
        ,   B_ {
                {0.3724,    0.9516,    0.2691},
                {0.1981,    0.9203,    0.4228},
                {0.4897,    0.0527,    0.5479},
                {0.3395,    0.7379,    0.9427}
            }
        ,   C_ {
                {0.1781,    0.9991,    0.0326,    0.8819},
                {0.1280,    0.1711,    0.5612,    0.6692}
            }
        ,   Q_ {
                {0.4447,    0.4192,    0.5073,    0.6125},
                {0.4192,    1.1026,    1.0150,    0.8676},
                {0.5073,    1.0150,    1.0083,    0.9772},
                {0.6125,    0.8676,    0.9772,    1.4596}
            }
        ,   R_ {
                {0.4442,    0.2368},
                {0.2368,    0.1547}
            }
        ,   xHat0_ {
                0.9436,    0.6377,    0.9577,    0.2407
            }
        ,   P0_ {
                {1.6604,    1.2726,    1.1100,    0.5482},
                {1.2726,    1.1217,    0.9461,    0.4314},
                {1.1100,    0.9461,    1.1303,    0.6500},
                {0.5482,    0.4314,    0.6500,    0.4815}
            }
        {
            kalman_.processNoiseCovariance(Q_);
            kalman_.measurementNoiseCovariance(R_);
            kalman_.stateEstimate(xHat0_);
            kalman_.stateCovariance(P0_);
        }


        static size_t constexpr NX = 4;
        static size_t constexpr NU = 3;
        static size_t constexpr NY = 2;

        KalmanFilter<Real> kalman_;
        Vec x0_;
        Vec u0_;
        Mat A_;
        Mat B_;
        Mat C_;
        Mat Q_;
        Mat R_;
        Vec xHat0_;
        Mat P0_;
    };


	TEST_F(KalmanFilterInitTest, testNx)
    {
        EXPECT_EQ(kalman_.nx(), NX);
    }


    TEST_F(KalmanFilterInitTest, testNy)
    {
        EXPECT_EQ(kalman_.ny(), NY);
    }


    TEST_F(KalmanFilterInitTest, testStateCovariance)
    {
        EXPECT_EQ(kalman_.stateCovariance(), blaze::ZeroMatrix<Real>(NX, NX));
    }


    TEST_F(KalmanFilterInitTest, testStateEstimate)
    {
        TMPC_EXPECT_EQ(kalman_.stateEstimate(), (blaze::ZeroVector<Real, blaze::columnVector>(NX)));
    }


    TEST_F(KalmanFilterInitTest, testProcessNoiseCovariance)
    {
        EXPECT_EQ(kalman_.processNoiseCovariance(), (blaze::ZeroMatrix<Real>(NX, NX)));
    }


    TEST_F(KalmanFilterInitTest, testMeasurementNoiseCovariance)
    {
        EXPECT_EQ(kalman_.measurementNoiseCovariance(), (blaze::ZeroMatrix<Real>(NY, NY)));
    }

    
    TEST_F(KalmanFilterTest, testUpdate)
    {
        Vec const y0 = C_ * x0_ + Vec {-0.2979, -0.1403} - C_ * kalman_.stateEstimate();
        kalman_.update(y0, C_);

        TMPC_EXPECT_APPROX_EQ(kalman_.stateEstimate(), (Vec {1.0672, 0.7747, 0.7956, 0.0509}), 1e-4, 0.);
        TMPC_EXPECT_APPROX_EQ(kalman_.stateCovariance(), (Mat {
            {0.4024,    0.2068,    0.0987,    0.0031},
            {0.2068,    0.2157,    0.1147,   -0.0068},
            {0.0987,    0.1147,    0.1118,    0.0205},
            {0.0031,   -0.0068,    0.0205,    0.0672}
        }), 1e-4, 0.);
    }


    TEST_F(KalmanFilterTest, testPredict1)
    {
        kalman_.predict(A_ * kalman_.stateEstimate() + B_ * u0_, A_);

        TMPC_EXPECT_APPROX_EQ(kalman_.stateEstimate(), (Vec {1.3990, 2.0599, 2.2287, 2.5951}), 1e-4, 0.);
        TMPC_EXPECT_APPROX_EQ(kalman_.stateCovariance(), (Mat {
            {1.3584,    2.3924,    2.2925,    2.3721},
            {2.3924,    5.4282,    4.8832,    4.6771},
            {2.2925,    4.8832,    4.8142,    4.6958},
            {2.3721,    4.6771,    4.6958,    5.0963}
        }), 1e-4, 0.);
    }


    TEST_F(KalmanFilterTest, testPredict2)
    {
        kalman_.predict(A_, B_, u0_);

        TMPC_EXPECT_APPROX_EQ(kalman_.stateEstimate(), (Vec {1.3990, 2.0599, 2.2287, 2.5951}), 1e-4, 0.);
        TMPC_EXPECT_APPROX_EQ(kalman_.stateCovariance(), (Mat {
            {1.3584,    2.3924,    2.2925,    2.3721},
            {2.3924,    5.4282,    4.8832,    4.6771},
            {2.2925,    4.8832,    4.8142,    4.6958},
            {2.3721,    4.6771,    4.6958,    5.0963}
        }), 1e-4, 0.);
    }


    TEST(KalmanFilterDifficultTest, testUpdate)
    {
        static size_t constexpr NX = 4;
        static size_t constexpr NW = 1;
        static size_t constexpr NH = 2;

                
        std::vector<StaticVector<NH>> const y_rec {
            StaticVector<NH> {-1.7889e-05, -1.30424e-05},
            StaticVector<NH> {-0.00022692, -0.000167837},
            StaticVector<NH> {-0.00316451, 0.000157017},
            StaticVector<NH> {-7.03035e-05, -1.82798e-05},
            StaticVector<NH> {-0.000253626, 0.000172559}
        };

        std::vector<StaticVector<NX + NW>> const x_rec {
            StaticVector<NX + NW> {0., 0., 0., 0.},
            StaticVector<NX + NW> {-1.29957e-05, 5.32864e-05, 1.24782e-05, -0.000849879, -5.21874e-10},
            StaticVector<NX + NW> {-0.000152557, -0.040238, 0.0001302, 0.00746183, -0.0063417},
            StaticVector<NX + NW> {0.000104791, 0.05479, -2.98354e-05, -0.00719635, 0.0179322},
            StaticVector<NX + NW> {-2.85142e-05, 0.0139012, 5.00306e-05, -0.00535151, 0.0307889}
        };

        std::vector<StaticMatrix<NX + NW, NX + NW>> const A_rec {
            StaticMatrix<NX + NW, NX + NW> {
                {   0.999998,  0.000997127,   0.00157084,  1.91313e-05,  3.98934e-07},
                {-0.00317067,     0.994295,      2.75675,    0.0342017,  0.000792177},
                {1.28311e-06,  2.30861e-06,     0.975111,  0.000696879, -3.20579e-07},
                {  0.0022518,   0.00405201,     -43.7332,     0.457453, -0.000562602},
                {          0,            0,            0,            0,            1},
            },
            StaticMatrix<NX + NW, NX + NW> {
                {   0.999998,  0.000997127,   0.00157084,  1.91313e-05,  3.98934e-07},
                {-0.00317067,     0.994295,      2.75675,    0.0342017,  0.000792177},
                {1.28311e-06,  2.30861e-06,     0.975111,  0.000696879, -3.20579e-07},
                {  0.0022518,   0.00405201,     -43.7332,     0.457453, -0.000562602},
                {          0,            0,            0,            0,            1},
            },
            StaticMatrix<NX + NW, NX + NW> {
                {   0.999998,  0.000997127,   0.00157084,  1.91313e-05,  3.98934e-07},
                {-0.00317067,     0.994295,      2.75675,    0.0342017,  0.000792177},
                {1.28311e-06,  2.30861e-06,     0.975111,  0.000696879, -3.20579e-07},
                {  0.0022518,   0.00405201,     -43.7332,     0.457453, -0.000562602},
                {          0,            0,            0,            0,            1},
            },
            StaticMatrix<NX + NW, NX + NW> {
                {   0.999998,  0.000997127,   0.00157084,  1.91313e-05,  3.98934e-07},
                {-0.00317067,     0.994295,      2.75675,    0.0342017,  0.000792177},
                {1.28311e-06,  2.30861e-06,     0.975111,  0.000696879, -3.20579e-07},
                {  0.0022518,   0.00405201,     -43.7332,     0.457453, -0.000562602},
                {          0,            0,            0,            0,            1},
            },
            StaticMatrix<NX + NW, NX + NW> {
                {   0.999998,  0.000997127,   0.00157084,  1.91313e-05,  3.98934e-07},
                {-0.00317067,     0.994295,      2.75675,    0.0342017,  0.000792177},
                {1.28311e-06,  2.30861e-06,     0.975111,  0.000696879, -3.20579e-07},
                {  0.0022518,   0.00405201,     -43.7332,     0.457453, -0.000562602},
                {          0,            0,            0,            0,            1},
            }
        };

        StaticMatrix<NH, NX + NW> C {
            {-100,            0,         -100,            0,            0},
            {   1,            0,            0,            0,            0}
        };
        
        // Create and prepare estimator.
        // EkfEstimator estimator(timeStepInSeconds);
        tmpc::KalmanFilter<Real> kalman(NX + NW, NH);
        kalman.stateEstimate(StaticVector<NX + NW> {0., 0., 0., 0., 0.});
        kalman.stateCovariance(StaticMatrix<NX + NW, NX + NW> {
            {1,            0,            0,            0,            0},
            {0,            1,            0,            0,            0},
            {0,            0,        10000,            0,            0},
            {0,            0,            0,        10000,            0},
            {0,            0,            0,            0,          100},
        });
        kalman.processNoiseCovariance(StaticMatrix<NX + NW, NX + NW> {
            {0,            0,            0,            0,            0},
            {0,            0,            0,            0,            0},
            {0,            0,            0,            0,            0},
            {0,            0,            0,            0,            0},
            {0,            0,            0,            0,         0.01},
        });
        kalman.measurementNoiseCovariance(StaticMatrix<NH, NH> {
            {1e-06,            0},
            {    0,        1e-06},
        });

        for (size_t i = 0; i < size(y_rec); ++i)
        {
            // x = integrator.integrate(model, x, u, w, timeStepInSeconds);
            StaticVector<NH> const y = y_rec[i];

            // Feed current control input to the estimator
            // DEBUG_predict(kalman, x_rec[i], A_rec[i]);
            kalman.predict(x_rec[i], A_rec[i]);

            // Update the estimate given the measurement
            // DEBUG_update(kalman, y_rec[i], C);
            kalman.update(y_rec[i], C);
        }
    }
}