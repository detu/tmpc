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


    // Test assignment behavior of symmetric matrices which are not exactly symmetric.
    TEST(SymmetricMatrixTest, testNotExactlySymmetricAssign)
    {
        blaze::SymmetricMatrix<blaze::DynamicMatrix<double>> const A {
            {1., 1e-9},
            {1e-10, 2.}
        };

        blaze::StaticMatrix<double, 2, 2, blaze::columnMajor> const B = A;
        // blaze::SymmetricMatrix<blaze::DynamicMatrix<double>> const B = A;
        TMPC_EXPECT_EQ(B, A);

        blaze::SymmetricMatrix<blaze::DynamicMatrix<double>> C {2ul};
        submatrix(C, 0, 0, 2, 2) = B;

        TMPC_EXPECT_EQ(C, B);
        TMPC_EXPECT_EQ(C, A);

        // std::clog << "A = \n" << A;
        // std::clog << "B = \n" << B;
        // std::clog << "C = \n" << C;
    }


    // Test init behavior of symmetric matrices which are not exactly symmetric.
    TEST(SymmetricMatrixTest, testNotExactlySymmetricInit)
    {
        blaze::SymmetricMatrix<blaze::DynamicMatrix<double>> const A {
            {1., 1e-9},
            {1e-10, 1.}
        };

        ASSERT_EQ(A(0, 1), 1e-9);
        ASSERT_EQ(A(1, 0), 1e-10);
    }


    TEST(MatrixTest, testInequality)
    {
        blaze::DynamicMatrix<double> A {
            {8.109899998657424e-02,     7.431769777268771e-02},
            {8.280816621049328e-02,     8.208946675009486e-03},
            {4.852287590951531e-02,     3.065557584184504e-02},
            {7.818775312987819e-02,     9.891525437807851e-03},
            {2.968117707753279e-02,     1.433470474925235e-02},
        };

        blaze::DynamicMatrix<double> B {
            {8.109899998657424e-02,     7.431769777368771e-02},
            {8.280816621049328e-02,     8.208946675009486e-03},
            {4.852287590951531e-02,     3.065557584184504e-02},
            {7.818775312987819e-02,     9.891525437807851e-03},
            {2.968117707753279e-02,     1.433470474925235e-02},
        };

        EXPECT_FALSE(exactEqual(A, B));
    }
}
