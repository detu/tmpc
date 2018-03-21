#include <tmpc/matrix/BlockMatrixView.hpp>
#include <tmpc/EigenKernel.hpp>

#include <gtest/gtest.h>


namespace tmpc :: testing
{
    class BlockMatrixViewTest 
    :   public ::testing::Test 
    {
    protected:
        using Kernel = EigenKernel<double>;
    };


    TEST_F(BlockMatrixViewTest, testCtor)
    {
        using Matrix = DynamicMatrix<Kernel>;
        
        Matrix m(10, 12);
        randomize(m);

        BlockMatrixView<Kernel, Matrix> block_view {m, {1, 3, 6}, {2, 7, 1, 2}};

        EXPECT_EQ(rows(block_view), 3);
        EXPECT_EQ(columns(block_view), 4);
    }

    
    TEST_F(BlockMatrixViewTest, testCtorInvalidColSize)
    {
        using Matrix = DynamicMatrix<Kernel>;
        
        Matrix m(10, 12);
        randomize(m);

        EXPECT_THROW((BlockMatrixView<Kernel, Matrix>(m, {1, 3, 6}, {2, 7, 1, 3})), std::invalid_argument);
    }


    TEST_F(BlockMatrixViewTest, testCtorInvalidRowSize)
    {
        using Matrix = DynamicMatrix<Kernel>;
        
        Matrix m(10, 12);
        randomize(m);

        EXPECT_THROW((BlockMatrixView<Kernel, Matrix>(m, {1, 2, 6}, {2, 7, 1, 2})), std::invalid_argument);
    }


    TEST_F(BlockMatrixViewTest, testBlockAccess)
    {
        using Matrix = DynamicMatrix<Kernel>;
        
        Matrix m(10, 12);
        randomize(m);

        std::array<size_t, 3> si = {1, 3, 6};
        std::array<size_t, 4> sj = {2, 7, 1, 2};
        BlockMatrixView<Kernel, Matrix> block_view {m, si, sj};

        size_t ii = 0;
        for (size_t i = 0; i < si.size(); ii += si[i++])
        {
            size_t jj = 0;
            for (size_t j = 0; j < sj.size(); jj += sj[j++])
                EXPECT_EQ(block_view(i, j), submatrix(m, ii, jj, si[i], sj[j]));
        }
    }
}