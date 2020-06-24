#include <tmpc/util/Unwrap.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    class UnwrapTest
    :   public Test
    {
    protected:
        Unwrap<double> unwrap_ {2.};
    };


    TEST_F(UnwrapTest, testFirstCall)
    {
        double val = 42.1;
        EXPECT_EQ(unwrap_(val), val);
    }


    TEST_F(UnwrapTest, testZeroDiff)
    {
        double val = 42.1;
        EXPECT_EQ(unwrap_(val), val);
        EXPECT_EQ(unwrap_(val), val);
    }


    TEST_F(UnwrapTest, testPositiveDiffNoWrap)
    {
        double v0 = 42.1;
        double v1 = v0 + 0.99 * unwrap_.halfPeriod();

        unwrap_(v0);
        EXPECT_EQ(unwrap_(v1), v1);
    }


    TEST_F(UnwrapTest, testPositiveDiffPositiveWrap)
    {
        double v0 = 42.1;
        double v1 = v0 + 0.99 * unwrap_.halfPeriod() + unwrap_.period();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 + 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST_F(UnwrapTest, testPositiveDiffPositiveWrap2)
    {
        double v0 = 42.1;
        double v1 = v0 + 0.99 * unwrap_.halfPeriod() + 2. * unwrap_.period();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 + 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST_F(UnwrapTest, testPositiveDiffNegativeWrap)
    {
        double v0 = 42.1;
        double v1 = v0 + 1.01 * unwrap_.halfPeriod();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 - 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST_F(UnwrapTest, testPositiveDiffNegativeWrap1)
    {
        double v0 = 42.1;
        double v1 = v0 + 1.01 * unwrap_.halfPeriod() + unwrap_.period();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 - 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST_F(UnwrapTest, testPositiveDiffNegativeWrap2)
    {
        double v0 = 42.1;
        double v1 = v0 + 1.01 * unwrap_.halfPeriod() + 2. * unwrap_.period();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 - 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST_F(UnwrapTest, testNegativeDiffNoWrap)
    {
        double v0 = 42.1;
        double v1 = v0 - 0.99 * unwrap_.halfPeriod();

        unwrap_(v0);
        EXPECT_EQ(unwrap_(v1), v1);
    }


    TEST_F(UnwrapTest, testNegativeDiffNegativeWrap)
    {
        double v0 = 0;
        double v1 = v0 - 0.99 * unwrap_.halfPeriod() - unwrap_.period();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 - 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST_F(UnwrapTest, testNegativeDiffNegativeWrap2)
    {
        double v0 = 42.1;
        double v1 = v0 - 0.99 * unwrap_.halfPeriod() - 2. * unwrap_.period();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 - 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST_F(UnwrapTest, testNegativeDiffPositiveWrap)
    {
        double v0 = 42.1;
        double v1 = v0 - 1.01 * unwrap_.halfPeriod();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 + 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST_F(UnwrapTest, testNegativeDiffPositiveWrap1)
    {
        double v0 = 42.1;
        double v1 = v0 - 1.01 * unwrap_.halfPeriod() - unwrap_.period();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 + 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST_F(UnwrapTest, testNegativeDiffPositiveWrap2)
    {
        double v0 = 42.1;
        double v1 = v0 - 1.01 * unwrap_.halfPeriod() - 2. * unwrap_.period();

        unwrap_(v0);
        EXPECT_NEAR(unwrap_(v1), v0 + 0.99 * unwrap_.halfPeriod(), 1e-10);
    }


    TEST(Unwrap2PiTest, test1)
    {
        Unwrap<double> unwrap(2. * M_PI);

        double v0 = 6.262908936000000;
        double v1 = 0.006951916032000;

        unwrap(v0);
        EXPECT_NEAR(unwrap(v1), 6.290137223211587, 1e-10);
    }
}