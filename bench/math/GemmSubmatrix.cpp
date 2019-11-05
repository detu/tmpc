#include <tmpc/Math.hpp>

#include <benchmark/benchmark.h>

#include <vector>


namespace tmpc :: benchmark
{
    using namespace ::benchmark;


    template <typename Real, size_t M>
    static void BM_gemm_static_submatrix(::benchmark::State& state)
    {
        size_t constexpr N = M;
        size_t constexpr K = M;
        
        blaze::StaticMatrix<Real, K, M, blaze::columnMajor> A;
        randomize(A);

        blaze::StaticMatrix<Real, K, N, blaze::columnMajor> B;
        randomize(B);

        blaze::StaticMatrix<Real, M, N, blaze::columnMajor> C;
        randomize(C);
        
        for (auto _ : state)
        {
            submatrix(C, 0, 0, M, N) += trans(submatrix(A, 0, 0, K, M)) * submatrix(B, 0, 0, K, N);
            // blaze::submatrix<0, 0, M, N>(C) += trans(blaze::submatrix<0, 0, K, M>(A)) * blaze::submatrix<0, 0, K, N>(B);
            ::benchmark::DoNotOptimize(A);
            ::benchmark::DoNotOptimize(B);
            ::benchmark::DoNotOptimize(C);
        }

        state.counters["flops"] = Counter(2 * M * N * K, Counter::kIsIterationInvariantRate);
        state.counters["m"] = M;
        state.counters["n"] = N;
        state.counters["k"] = K;
    }


    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 1);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 2);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 3);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 4);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 5);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 6);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 7);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 8);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 9);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 10);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 11);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 12);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 13);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 14);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 15);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 16);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 17);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 18);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 19);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 20);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 21);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 22);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 23);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 24);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 25);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 26);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 27);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 28);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 29);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 30);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 31);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 32);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 33);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 34);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 35);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 36);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 37);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 38);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 39);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 40);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 41);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 42);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 43);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 44);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 45);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 46);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 47);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 48);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 49);
    BENCHMARK_TEMPLATE(BM_gemm_static_submatrix,double, 50);
}
