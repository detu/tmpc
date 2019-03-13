#include <tmpc/estimation/KalmanFilter.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    using Real = double;
    using Mat = blaze::DynamicMatrix<Real>;
    using Vec = blaze::DynamicVector<Real, blaze::columnVector>;


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
}