#include <tmpc/Matrix.hpp>

#include <gtest/gtest.h>

using namespace tmpc;

template <typename MT, bool SO>
std::istream& operator>>(std::istream& is, Matrix<MT, SO>& m)
{
	for (size_t i = 0; i < rows(m); ++i)
		for (size_t j = 0; j < columns(m); ++j)
			; //is >> (~m)(i, j);

	return is;
}

template <typename VT, bool TF>
std::istream& operator>>(std::istream& is, Vector<VT, TF>& v)
{
	for (size_t i = 0; i < size(v); ++i)
		; //is >> (~v)[i];

	return is;
}

/*
TEST(TestMatrix, InputStaticMatrixTest)
{
    typedef StaticMatrix<double, 3, 4> M;
    M m;
    std::cin >> static_cast<Matrix<M, rowMajor>&>(m);
}
*/

TEST(TestMatrix, InputStaticVectorTest)
{
    typedef StaticVector<double, 3> V;
    V v;
    std::cin >> static_cast<Vector<V, columnVector>&>(m);
}

#ifndef USE_BLAZE

TEST(TestMatrix, IsVectorTest)
{
    static_assert(IsVector<StaticVector<double, 2>>::value, 
        "IsVector<StaticVector>::value must be true!");
}

TEST(TestMatrix, SizeTest)
{
    static_assert(Size<StaticVector<double, 2, false>>::value == 2, 
        "Size<StaticVector<double, 2>>::value must be equal to 2!");
}

#endif