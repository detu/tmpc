#include <tmpc/Math.hpp>

#include <benchmark/benchmark.h>

#include <vector>


namespace tmpc :: benchmark
{
    using namespace ::benchmark;


    template <typename Real, size_t M, size_t N, size_t K>
    static void BM_gemm_blaze_static(::benchmark::State& state)
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

        state.counters["flops"] = Counter(M * N * K, Counter::kIsIterationInvariantRate);
        state.counters["m"] = M;
        state.counters["n"] = N;
        state.counters["k"] = K;
    }


    template <typename Real>
    static void BM_gemm_blaze_dynamic(::benchmark::State& state)
    {
        size_t const m = state.range(0);

        blaze::DynamicMatrix<Real, blaze::columnMajor> A(m, m);
        randomize(A);

        blaze::DynamicMatrix<Real, blaze::columnMajor> B(m, m);
        randomize(B);

        blaze::DynamicMatrix<Real, blaze::columnMajor> C(m, m);
        randomize(C);

        blaze::DynamicMatrix<Real, blaze::columnMajor> D(m, m);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(D = C + trans(A) * B);

        state.counters["flops"] = Counter(m * m * m, Counter::kIsIterationInvariantRate);
        state.counters["m"] = m;
    }


    template <typename Real>
    static void BM_gemm_loop_naive(::benchmark::State& state)
    {
        size_t const m = state.range(0);

        std::vector<Real> A(m * m), B(m * m), C(m * m), D(m * m);
        double * pA = A.data(), * pB = B.data(), * pC = C.data(), * pD = D.data();
        
        for (auto _ : state)
        {
            for (size_t j = 0; j < m; ++j)
                for (size_t i = 0; i < m; ++i)
                {
                    double s = pC[i + j * m];

                    for (size_t k = 0; k < m; ++k)
                        s += pA[k + i * m] * pB[k + j * m];

                    pD[i + j * m] = s;
                }
        }

        state.counters["flops"] = Counter(m * m * m, Counter::kIsIterationInvariantRate);
        state.counters["m"] = m;
    }


    template <typename Real>
    static void BM_gemm_loop_optimized(::benchmark::State& state)
    {
        size_t const m = state.range(0);

        std::vector<Real> A(m * m), B(m * m), C(m * m), D(m * m);
        double const * pA = A.data();
        double const * pB = B.data();
        double const * pC = C.data();
        double * pD = D.data();
        
        for (auto _ : state)
        {
            for (size_t j = 0; j < m; ++j)
                for (size_t i = 0; i < m; ++i)
                {
                    double s = pC[i + j * m];

                    double const * pa = pA + i * m;
                    double const * pb = pB + j * m;

                    for (size_t k = 0; k < m; ++k)
                        s += pa[k] * pb[k];

                    pD[i + j * m] = s;
                }
        }

        state.counters["flops"] = Counter(m * m * m, Counter::kIsIterationInvariantRate);
        state.counters["m"] = m;
    }


    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 2, 2, 2);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 3, 3, 3);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 5, 5, 5);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 10, 10, 10);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 20, 20, 20);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 30, 30, 30);

    BENCHMARK_TEMPLATE(BM_gemm_blaze_dynamic, double)->DenseRange(1, 40);
    
    BENCHMARK_TEMPLATE(BM_gemm_loop_naive, double)->DenseRange(1, 40);
    BENCHMARK_TEMPLATE(BM_gemm_loop_optimized, double)->DenseRange(1, 40);
}
