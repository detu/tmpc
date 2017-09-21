#include <tmpc/EigenKernel.hpp>

#include <gtest/gtest.h>

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

namespace tmpc :: eigen_adaptor :: testing
{
    using Kernel = EigenKernel<double>;

    TEST(TestMatrix, testInputStaticMatrix)
    {
        typedef StaticMatrix<Kernel, 2, 3> M;
        M m;

        std::istringstream is("1 2 3\n 4 5 6");
        is >> m;

        EXPECT_EQ(m, (M {{1., 2., 3.}, {4., 5., 6.}}));
    }

    TEST(TestMatrix, testInputStaticVector)
    {
        typedef StaticVector<Kernel, 3> V;
        V v;

        std::istringstream is("1 2 3");
        is >> v;

        EXPECT_EQ(v, (V {1., 2., 3.}));
    }

    TEST(TestMatrix, testIsVector)
    {
        static_assert(IsVector<StaticVector<Kernel, 2>>::value, 
            "IsVector<StaticVector>::value must be true!");
    }

    TEST(TestMatrix, testSize)
    {
        static_assert(Size<StaticVector<Kernel, 2, columnVector>>::value == 2, 
            "Size<StaticVector<double, 2>>::value must be equal to 2!");
    }

    TEST(TestMatrix, testColumns)
    {
        DynamicMatrix<Kernel> m(10, 20);
        EXPECT_EQ(columns(m), 20);
    }

    TEST(TestMatrix, testRows)
    {
        DynamicMatrix<Kernel> m(10, 20);
        EXPECT_EQ(rows(m), 10);
    }

    TEST(TestMatrix, testIsDerivedFromMatrixTag)
    {
        static_assert(std::is_base_of<MatrixTag, DynamicMatrix<Kernel>>::value, 
            "DynamicMatrix must be derived from MatrixTag");
    }
}