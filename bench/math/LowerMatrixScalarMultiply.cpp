#include <blaze/Math.h>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Real, size_t N, bool SO>
    static Real LowerMatrixScalarMultiplyStatic_Impl(blaze::StaticMatrix<Real, N, N, SO>& A);


    template <typename Real, size_t N, bool SO>
    static Real LowerMatrixScalarMultiplyStaticLoop_Impl(blaze::StaticMatrix<Real, N, N, SO>& A);


    template <typename Real, size_t N, bool SO>
    static void BM_LowerMatrixScalarMultiplyStatic(::benchmark::State& state)
    {
        blaze::StaticMatrix<Real, N, N, SO> A;        
        randomize(A);

        for (auto _ : state)
            ::benchmark::DoNotOptimize(LowerMatrixScalarMultiplyStatic_Impl(A));
    }


    template <typename Real, size_t N, bool SO>
    static void BM_LowerMatrixScalarMultiplyStaticLoop(::benchmark::State& state)
    {
        blaze::StaticMatrix<Real, N, N, SO> A;        
        randomize(A);

        for (auto _ : state)
            ::benchmark::DoNotOptimize(LowerMatrixScalarMultiplyStaticLoop_Impl(A));
    }


    template <typename Real, size_t N, bool SO>
    static Real LowerMatrixScalarMultiplyStatic_Impl(blaze::StaticMatrix<Real, N, N, SO>& A)
    {
        for (size_t k = 0; k < N; ++k)
        {
            size_t const rs = N - k - 1;
            auto A21 = submatrix(A, k + 1, k, rs, 1);

            A21 *= 1.1;
        }

        return A(N - 1, N - 1);

        // auto A21 = submatrix(A, 1, 0, N - 1, 1);
        // A21 *= 1.1;

        // return A21(N - 2, 0);
    }


    template <typename Real, size_t N, bool SO>
    static Real LowerMatrixScalarMultiplyStaticLoop_Impl(blaze::StaticMatrix<Real, N, N, SO>& A)
    {
        for (size_t k = 0; k < N; ++k)
        {
            size_t const rs = N - k - 1;
            auto A21 = submatrix(A, k + 1, k, rs, 1);

            for (size_t i = 0; i < rs; ++i)
                A21(i, 0) *= 1.1;
        }

        return A(N - 1, N - 1);

        // auto A21 = submatrix(A, 1, 0, N - 1, 1);

        // for (size_t i = 0; i < N - 1; ++i)
        //     A21(i, 0) *= 1.1;

        // return A(N - 2, 0);
    }
    
    
    BENCHMARK_TEMPLATE(BM_LowerMatrixScalarMultiplyStatic, double, 5, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_LowerMatrixScalarMultiplyStaticLoop, double, 5, blaze::columnMajor);
}
