#include <tmpc/Matrix.hpp>

#include <gtest/gtest.h>

using namespace tmpc;

template <typename MT, bool SO>
std::istream& operator>>(std::istream& is, Matrix<MT, SO>& m)
{
	for (size_t i = 0; i < rows(m); ++i)
		for (size_t j = 0; j < columns(m); ++j)
			is >> (~m)(i, j);

	return is;
}

template <typename VT, bool TF>
std::istream& operator>>(std::istream& is, Vector<VT, TF>& v)
{
	for (size_t i = 0; i < size(v); ++i)
		is >> (~v)[i];

	return is;
}

TEST(TestMatrix, testInputStaticMatrix)
{
    typedef StaticMatrix<double, 2, 3> M;
    M m;

    std::istringstream is("1 2 3\n 4 5 6");
    is >> m;

    EXPECT_EQ(m, (M {{1., 2., 3.}, {4., 5., 6.}}));
}

TEST(TestMatrix, testInputStaticVector)
{
    typedef StaticVector<double, 3> V;
    V v;

    std::istringstream is("1 2 3");
    is >> v;

    EXPECT_EQ(v, (V {1., 2., 3.}));
}

#ifndef USE_BLAZE

TEST(TestMatrix, testIsVector)
{
    static_assert(IsVector<StaticVector<double, 2>>::value, 
        "IsVector<StaticVector>::value must be true!");
}

TEST(TestMatrix, testSize)
{
    static_assert(Size<StaticVector<double, 2, false>>::value == 2, 
        "Size<StaticVector<double, 2>>::value must be equal to 2!");
}

TEST(TestMatrix, testDynamicSubmatrixConstructorFromEigenBlock)
{
    using EigenMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    DynamicMatrix<double> m(20, 30);
    Eigen::Block<EigenMatrix, Eigen::Dynamic, Eigen::Dynamic> b(m, 2, 3, 4, 5);

    Submatrix<DynamicMatrix<double>> sm(b);
}

#endif