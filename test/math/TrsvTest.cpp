#include <tmpc/math/Trsv.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    TEST(TrsvTest, testLowerRowMajor)
    {
        size_t constexpr N = 3;
        
        blaze::LowerMatrix<blaze::StaticMatrix<double, N, N, blaze::rowMajor>> A;
        randomize(A);

        blaze::StaticVector<double, N, blaze::columnVector> b, x;
        randomize(b);

        tmpc::trsv(A, b, x);

        TMPC_EXPECT_APPROX_EQ(x, evaluate(inv(A) * b), 1e-10, 1e-10);
    }


    TEST(TrsvTest, testLowerColumnMajor)
    {
        size_t constexpr N = 3;
        
        blaze::LowerMatrix<blaze::StaticMatrix<double, N, N, blaze::columnMajor>> A;
        randomize(A);

        blaze::StaticVector<double, N, blaze::columnVector> b, x;
        randomize(b);

        tmpc::trsv(A, b, x);

        TMPC_EXPECT_APPROX_EQ(x, evaluate(inv(A) * b), 1e-10, 1e-10);
    }


    TEST(TrsvTest, testUpperRowMajor)
    {
        size_t constexpr N = 3;
        
        blaze::UpperMatrix<blaze::StaticMatrix<double, N, N, blaze::rowMajor>> A;
        randomize(A);

        blaze::StaticVector<double, N, blaze::columnVector> b, x;
        randomize(b);

        tmpc::trsv(A, b, x);

        TMPC_EXPECT_APPROX_EQ(x, evaluate(inv(A) * b), 1e-10, 1e-10);
    }


    TEST(TrsvTest, testUpperColumnMajor)
    {
        size_t constexpr N = 3;
        
        blaze::UpperMatrix<blaze::StaticMatrix<double, N, N, blaze::columnMajor>> A;
        randomize(A);

        blaze::StaticVector<double, N, blaze::columnVector> b, x;
        randomize(b);

        tmpc::trsv(A, b, x);

        TMPC_EXPECT_APPROX_EQ(x, evaluate(inv(A) * b), 1e-10, 1e-10);
    }


    TEST(TrsvTest, testLowerRowMajorTrans)
    {
        size_t constexpr N = 3;
        
        blaze::LowerMatrix<blaze::StaticMatrix<double, N, N, blaze::rowMajor>> A;
        randomize(A);

        blaze::StaticVector<double, N, blaze::rowVector> b, x;
        randomize(b);

        tmpc::trsv(b, A, x);

        TMPC_EXPECT_APPROX_EQ(x, evaluate(b * inv(A)), 1e-10, 1e-10);
    }


    TEST(TrsvTest, testLowerColumnMajorTrans)
    {
        size_t constexpr N = 3;
        
        blaze::LowerMatrix<blaze::StaticMatrix<double, N, N, blaze::columnMajor>> A;
        randomize(A);

        blaze::StaticVector<double, N, blaze::rowVector> b, x;
        randomize(b);

        tmpc::trsv(b, A, x);

        TMPC_EXPECT_APPROX_EQ(x, evaluate(b * inv(A)), 1e-10, 1e-10);
    }


    TEST(TrsvTest, testUpperRowMajorTrans)
    {
        size_t constexpr N = 3;
        
        blaze::UpperMatrix<blaze::StaticMatrix<double, N, N, blaze::rowMajor>> A;
        randomize(A);

        blaze::StaticVector<double, N, blaze::rowVector> b, x;
        randomize(b);

        tmpc::trsv(b, A, x);

        TMPC_EXPECT_APPROX_EQ(x, evaluate(b * inv(A)), 1e-10, 1e-10);
    }


    TEST(TrsvTest, testUpperColumnMajorTrans)
    {
        size_t constexpr N = 3;
        
        blaze::UpperMatrix<blaze::StaticMatrix<double, N, N, blaze::columnMajor>> A;
        randomize(A);

        blaze::StaticVector<double, N, blaze::rowVector> b, x;
        randomize(b);

        tmpc::trsv(b, A, x);

        TMPC_EXPECT_APPROX_EQ(x, evaluate(b * inv(A)), 1e-10, 1e-10);
    }
}
