#include <tmpc/Matrix.hpp>

#include <gtest/gtest.h>

namespace A
{
    template <typename T>
    struct AA
    {

    };

    template <typename T>
    struct BB
    :   AA<BB<T>>
    {        
    };

    using namespace Eigen;
}

namespace A
{
    template <typename T>
    void f(AA<T> const&)
    {
        std::cout << "OK!" << std::endl;
    }
}

namespace Eigen
{
    template <typename T>
    void g(MatrixBase<T> const&)
    {
        std::cout << "OK!" << std::endl;
    }

    template <typename T>
    void lpNorm(MatrixBase<T> const&, unsigned P)
    {
        std::cout << "OK!" << std::endl;
    }
}

TEST(TestMatrix, testLpNorm)
{
    A::BB<int> a;
    f(a);

    Eigen::MatrixXd m;
    g(m);

    lpNorm(m, 1);
}

using namespace tmpc;

template <typename MT, StorageOrder SO>
std::istream& operator>>(std::istream& is, Matrix<MT, SO>& m)
{
	for (size_t i = 0; i < rows(m); ++i)
		for (size_t j = 0; j < columns(m); ++j)
			is >> (~m)(i, j);

	return is;
}

template <typename VT, TransposeFlag TF>
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

TEST(TestMatrix, testColumns)
{
    DynamicMatrix<double> m(10, 20);
    EXPECT_EQ(columns(m), 20);
}

TEST(TestMatrix, testRows)
{
    DynamicMatrix<double> m(10, 20);
    EXPECT_EQ(rows(m), 10);
}

TEST(TestMatrix, testIsDerivedFromMatrixTag)
{
    static_assert(std::is_base_of<MatrixTag, DynamicMatrix<double>>::value, 
        "DynamicMatrix must be derived from MatrixTag");
}

#endif