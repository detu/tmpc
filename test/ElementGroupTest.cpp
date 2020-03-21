#include <tmpc/ElementGroup.hpp>
#include <tmpc/Testing.hpp>

#include <blaze/Math.h>


namespace tmpc :: testing
{
    TEST(ElementGroupTest, testDefineMultipleElementGroups)
    {
        auto constexpr group1 = ElementGroup<0, 2> {};
        auto constexpr group2 = ElementGroup<group1.end, 3> {};

        EXPECT_EQ(group2.end, 5);
    }


    TEST(ElementGroupTest, testVectorGet)
    {
        blaze::DynamicVector<double> x(5);
        randomize(x);

        ElementGroup<1, 2> constexpr comp;        
        EXPECT_TRUE(exactEqual(subvector(x, comp), subvector(x, 1, 2)));
    }


    TEST(ElementGroupTest, testVectorAssign)
    {
        blaze::DynamicVector<double> x(5);
        randomize(x);

        ElementGroup<1, 2> constexpr comp;
        blaze::DynamicVector<double> b(2);
        randomize(b);
        subvector(x, comp) = b;

        EXPECT_TRUE(exactEqual(subvector(x, 1, 2), b));
    }


    TEST(ElementGroupTest, testMatrixGet)
    {
        blaze::DynamicMatrix<double> m(5, 6);
        randomize(m);

        ElementGroup<1, 2> constexpr comp1;
        ElementGroup<3, 3> constexpr comp2;
        
        EXPECT_TRUE(exactEqual(submatrix(m, comp1, comp2), submatrix(m, 1, 3, 2, 3)));
    }


    TEST(ElementGroupTest, testMatrixAssign)
    {
        blaze::DynamicMatrix<double> m(5, 6);
        randomize(m);
        
        ElementGroup<1, 2> constexpr comp1;
        ElementGroup<3, 3> constexpr comp2;

        blaze::DynamicMatrix<double> b(2, 3);
        randomize(b);
        submatrix(m, comp1, comp2) = b;

        EXPECT_TRUE(exactEqual(submatrix(m, 1, 3, 2, 3), b));
    }
}
