#include <tmpc/Math.hpp>

#include <benchmark/benchmark.h>

#include <vector>


namespace tmpc :: benchmark
{
    using namespace ::benchmark;


    template <typename Real, size_t M>
    static void BM_gemm_blaze_static(::benchmark::State& state)
    {
        size_t constexpr N = M;
        size_t constexpr K = M;
        
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
    static void BM_gemm_blas(::benchmark::State& state)
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
            gemm(C, trans(A), B, 1.0, 1.0);

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


    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 1);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 2);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 3);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 4);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 5);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 6);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 7);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 8);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 9);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 10);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 11);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 12);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 13);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 14);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 15);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 16);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 17);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 18);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 19);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 20);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 21);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 22);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 23);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 24);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 25);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 26);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 27);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 28);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 29);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 30);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 31);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 32);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 33);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 34);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 35);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 36);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 37);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 38);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 39);
    BENCHMARK_TEMPLATE(BM_gemm_blaze_static, double, 40);    

    BENCHMARK_TEMPLATE(BM_gemm_blaze_dynamic, double)->DenseRange(1, 40);

    BENCHMARK_TEMPLATE(BM_gemm_blas, double)->DenseRange(1, 40);
    
    BENCHMARK_TEMPLATE(BM_gemm_loop_naive, double)->DenseRange(1, 40);
    BENCHMARK_TEMPLATE(BM_gemm_loop_optimized, double)->DenseRange(1, 40);
}
