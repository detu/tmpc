#include <tmpc/Matrix.hpp>

#include <gtest/gtest.h>

using namespace tmpc;

class DynamicSubmatrixTest : public ::testing::Test 
{
protected:
    using Matrix = DynamicMatrix<double, rowMajor>;

    Matrix m {20, 30};
};

TEST_F(DynamicSubmatrixTest, testConstructorFromEigenBlockRRef)
{
    Submatrix<Matrix> sm(m.block(2, 3, 4, 5));

    EXPECT_EQ(m.IsRowMajor, sm.IsRowMajor);
}

TEST_F(DynamicSubmatrixTest, testConstructorFromEigenBlockConstRef)
{
    auto const b = m.block(2, 3, 4, 5);

    Submatrix<Matrix> sm(b);

    EXPECT_EQ(m.IsRowMajor, b.IsRowMajor);
}

TEST_F(DynamicSubmatrixTest, testGet)
{
    using Matrix = DynamicMatrix<double, rowMajor>;

    Matrix m(20, 30);
    
    EXPECT_EQ(submatrix(m, 2, 3, 4, 5), m.block(2, 3, 4, 5));
}