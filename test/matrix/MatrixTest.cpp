#include <blaze/Math.h>

#include <tmpc/Testing.hpp>

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
    /// Note that depending on the storage order of the system matrix and the given right-hand side 
    /// the trsv() function solves different equation systems:
    ///     Single right-hand side:
    ///     A *x=b if A is column-major
    ///     A^T*x=b if A is row-major
    /// (see https://bitbucket.org/blaze-lib/blaze/wiki/LAPACK%20Functions#!linear-system-solver)
    TEST(MatrixTest, testTrsvRowMajor)
    {
        blaze::DynamicMatrix<double, blaze::rowMajor> Lambda {
            {1., std::numeric_limits<double>::signaling_NaN()},
            {-2., 3.}
        };

        blaze::DynamicVector<double> l {4., 7.};

        trsv(Lambda, l, 'L', 'T', 'N');

        blaze::DynamicVector<double> expected {4., 5.};

        EXPECT_EQ(forcePrint(l), forcePrint(expected));
    }
}
