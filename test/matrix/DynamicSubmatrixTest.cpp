#include <tmpc/Matrix.hpp>

#include <gtest/gtest.h>

using namespace tmpc;

class DynamicSubmatrixTest : public ::testing::Test 
{
protected:
    using Matrix = DynamicMatrix<double, rowMajor>;
};

TEST_F(DynamicSubmatrixTest, testConstructorFromEigenBlockRRef)
{
    Matrix m(20, 30);

    Submatrix<Matrix> sm(m.block(2, 3, 4, 5));

    EXPECT_EQ(m.IsRowMajor, sm.IsRowMajor);
}

TEST_F(DynamicSubmatrixTest, testConstructorFromEigenBlockConstRef)
{
    Matrix m(20, 30);
    auto const b = m.block(2, 3, 4, 5);

    Submatrix<Matrix> sm(b);
    
    EXPECT_EQ(m.IsRowMajor, b.IsRowMajor);
}

TEST_F(DynamicSubmatrixTest, testConstructorFromConstSubmatrix)
{
    Matrix const m(20, 30);
    Submatrix<Matrix const> sm(submatrix(m, 2, 3, 4, 5));
    EXPECT_EQ(sm, m.block(2, 3, 4, 5));
}

TEST_F(DynamicSubmatrixTest, testGet)
{
    Matrix m(20, 30);
    EXPECT_EQ(submatrix(m, 2, 3, 4, 5), m.block(2, 3, 4, 5));
}

TEST_F(DynamicSubmatrixTest, testGetConst)
{
    Matrix const m(20, 30);
    EXPECT_EQ(submatrix(m, 2, 3, 4, 5), m.block(2, 3, 4, 5));
}

TEST_F(DynamicSubmatrixTest, testScalarAssignment)
{
    Matrix m(4, 5, 0.);
    submatrix(m, 1, 2, 2, 3) = 1.;
    EXPECT_EQ(m, (Matrix {
        {0, 0, 0, 0, 0},
        {0, 0, 1, 1, 1},
        {0, 0, 1, 1, 1},
        {0, 0, 0, 0, 0}
    }));
}