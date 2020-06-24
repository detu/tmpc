#include <tmpc/math/Rank.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    TEST(RankTest, testFullRank)
    {
        blaze::StaticMatrix<double, 2, 3> m {
            {1., 2., 3.},
            {4., 5., 6.}
        };

        EXPECT_EQ(rank(m), 2);
    }


    TEST(RankTest, testDeficientRowRank)
    {
        blaze::StaticMatrix<double, 2, 3> m {
            {1., 2., 3.},
            {2., 4., 6.}
        };

        EXPECT_EQ(rank(m), 1);
    }


    TEST(RankTest, testDeficientColumnRank)
    {
        blaze::StaticMatrix<double, 3, 2> m {
            {1., 2.},
            {2., 4.},
            {3., 6.}
        };

        EXPECT_EQ(rank(m), 1);
    }
}