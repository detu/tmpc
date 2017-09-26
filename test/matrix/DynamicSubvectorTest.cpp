#include <tmpc/EigenKernel.hpp>

#include <gtest/gtest.h>

using namespace tmpc;

using Kernel = EigenKernel<double>;

class DynamicSubvectorTest : public ::testing::Test 
{
protected:
    using Vector = DynamicVector<Kernel>;
};

TEST_F(DynamicSubvectorTest, testConstructorFromEigenBlockRRef)
{
    Vector v(20);
    Subvector<Kernel, Vector> sv(v.segment(2, 4));

    EXPECT_EQ(v.IsRowMajor, sv.IsRowMajor);
}

TEST_F(DynamicSubvectorTest, testConstructorFromEigenBlockConstRef)
{
    Vector v(20);
    auto const b = v.segment(2, 4);

    Subvector<Kernel, Vector> sv(b);

    EXPECT_EQ(sv.IsRowMajor, b.IsRowMajor);
}

TEST_F(DynamicSubvectorTest, testConstructorFromConstSubvector)
{
    Vector const v(20);
    Subvector<Kernel, Vector const> sv(subvector(v, 2, 4));
    EXPECT_EQ(sv, v.segment(2, 4));
}

TEST_F(DynamicSubvectorTest, testConstructorFromSubvector)
{
    Vector v(20);
    Subvector<Kernel, Vector> sv(subvector(v, 2, 4));
    EXPECT_EQ(sv, v.segment(2, 4));
}

TEST_F(DynamicSubvectorTest, testGet)
{
    Vector v(20);
    EXPECT_EQ(subvector(v, 2, 4), v.segment(2, 4));
}

TEST_F(DynamicSubvectorTest, testGetConst)
{
    Vector const v(20);
    EXPECT_EQ(subvector(v, 2, 4), v.segment(2, 4));
}

TEST_F(DynamicSubvectorTest, testScalarAssignment)
{
    Vector v(20);
    subvector(v, 2, 4) = 42.;
    EXPECT_EQ(v[2], 42);
    EXPECT_EQ(v[3], 42);
    EXPECT_EQ(v[4], 42);
    EXPECT_EQ(v[5], 42);
}