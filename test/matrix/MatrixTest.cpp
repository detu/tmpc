#include <blaze/Math.h>

#include <tmpc/test_tools.hpp>

#include <limits>


namespace tmpc :: testing
{
    /// @brief Test trsv() for column-major matrices.
    TEST(MatrixTest, testTrsvColumnMajor)
    {
        blaze::DynamicMatrix<double, blaze::columnMajor> Lambda {
            {1., std::numeric_limits<double>::signaling_NaN()},
            {-2., 3.}
        };

        blaze::DynamicVector<double> l {4., 7.};
        
        trsv(Lambda, l, 'L', 'N', 'N');

        blaze::DynamicVector<double> expected {4., 5.};
        
        EXPECT_EQ(forcePrint(l), forcePrint(expected));
    }


    /// @brief Test trsv() for row-major matrices.
    ///
    /// Reprodces the following Blaze bug: https://bitbucket.org/blaze-lib/blaze/issues/216
    TEST(MatrixTest, testTrsvRowMajor)
    {
        blaze::DynamicMatrix<double, blaze::rowMajor> Lambda {
            {1., std::numeric_limits<double>::signaling_NaN()},
            {-2., 3.}
        };

        blaze::DynamicVector<double> l {4., 7.};

        trsv(Lambda, l, 'L', 'N', 'N');

        blaze::DynamicVector<double> expected {4., 5.};

        EXPECT_EQ(forcePrint(l), forcePrint(expected));
    }
}
