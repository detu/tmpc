#include <tmpc/estimation/KalmanFilter.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    using Real = double;


    class KalmanFilterInitTest
    :   public Test
    {
    protected:
        static size_t constexpr NX = 4;
        static size_t constexpr NU = 3;
        static size_t constexpr NY = 2;

        KalmanFilter<Real> kalman_ {NX, NU, NY};
    };


	TEST_F(KalmanFilterInitTest, testNx)
    {
        EXPECT_EQ(kalman_.nx(), NX);
    }


    TEST_F(KalmanFilterInitTest, testNu)
    {
        EXPECT_EQ(kalman_.nu(), NU);
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


    TEST_F(KalmanFilterInitTest, testA)
    {
        EXPECT_EQ(kalman_.A(), (blaze::IdentityMatrix<Real>(NX)));
    }


    TEST_F(KalmanFilterInitTest, testB)
    {
        EXPECT_EQ(kalman_.B(), (blaze::ZeroMatrix<Real>(NX, NU)));
    }


    TEST_F(KalmanFilterInitTest, testC)
    {
        TMPC_EXPECT_EQ(kalman_.C(), (blaze::ZeroMatrix<Real>(NY, NX)));
    }


    TEST_F(KalmanFilterInitTest, testProcessNoiseCovariance)
    {
        EXPECT_EQ(kalman_.processNoiseCovariance(), (blaze::ZeroMatrix<Real>(NX, NX)));
    }


    TEST_F(KalmanFilterInitTest, testMeasurementNoiseCovariance)
    {
        EXPECT_EQ(kalman_.measurementNoiseCovariance(), (blaze::ZeroMatrix<Real>(NY, NY)));
    }
}