#include <tmpc/Math.hpp>

#include <benchmark/benchmark.h>

#include <vector>


namespace tmpc :: benchmark
{
    template <typename Real, size_t M, size_t N, size_t K>
    static void BM_gemm_blaze(::benchmark::State& state)
    {
        blaze::StaticMatrix<Real, K, M, blaze::columnMajor> A;
        randomize(A);

        blaze::StaticMatrix<Real, K, N, blaze::columnMajor> B;
        randomize(B);

        blaze::StaticMatrix<Real, M, N, blaze::columnMajor> C;
        randomize(C);

        blaze::StaticMatrix<Real, M, N, blaze::columnMajor> D;
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(D = C + trans(A) * B);
    }


    BENCHMARK_TEMPLATE(BM_gemm_blaze, double, 2, 2, 2);
    BENCHMARK_TEMPLATE(BM_gemm_blaze, double, 3, 3, 3);
    BENCHMARK_TEMPLATE(BM_gemm_blaze, double, 5, 5, 5);
    BENCHMARK_TEMPLATE(BM_gemm_blaze, double, 10, 10, 10);
    BENCHMARK_TEMPLATE(BM_gemm_blaze, double, 20, 20, 20);
    BENCHMARK_TEMPLATE(BM_gemm_blaze, double, 30, 30, 30);
}
