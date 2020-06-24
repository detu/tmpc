#include <blaze/Math.h>
#include <Eigen/Dense>

#include <tmpc/math/Trsv.hpp>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Real, size_t M, bool SO>
    static void BM_trsv_Lower_tmpc_Static(::benchmark::State& state)
    {
        blaze::LowerMatrix<blaze::StaticMatrix<Real, M, M, SO>> A;        
        blaze::StaticVector<Real, M, blaze::columnVector> b;
        blaze::StaticVector<Real, M, blaze::columnVector> x;

        randomize(A);
        randomize(b);
        
        for (auto _ : state)
        {
            tmpc::trsv(A, b, x);
            ::benchmark::DoNotOptimize(x[M - 1]);
        }
    }


    template <typename Real, size_t M, bool SO>
    static void BM_trsv_TransLower_tmpc_Static(::benchmark::State& state)
    {
        blaze::LowerMatrix<blaze::StaticMatrix<Real, M, M, SO>> A;        
        blaze::StaticVector<Real, M, blaze::rowVector> b;
        blaze::StaticVector<Real, M, blaze::rowVector> x;

        randomize(A);
        randomize(b);
        
        for (auto _ : state)
        {
            tmpc::trsv(b, A, x);
            ::benchmark::DoNotOptimize(x[M - 1]);
        }
    }


    BENCHMARK_TEMPLATE(BM_trsv_Lower_tmpc_Static, double, 1, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_trsv_Lower_tmpc_Static, double, 4, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_trsv_Lower_tmpc_Static, double, 35, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_trsv_Lower_tmpc_Static, double, 60, blaze::rowMajor);

    BENCHMARK_TEMPLATE(BM_trsv_Lower_tmpc_Static, double, 1, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_trsv_Lower_tmpc_Static, double, 4, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_trsv_Lower_tmpc_Static, double, 35, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_trsv_Lower_tmpc_Static, double, 60, blaze::columnMajor);

    BENCHMARK_TEMPLATE(BM_trsv_TransLower_tmpc_Static, double, 1, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_trsv_TransLower_tmpc_Static, double, 4, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_trsv_TransLower_tmpc_Static, double, 35, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_trsv_TransLower_tmpc_Static, double, 60, blaze::rowMajor);

    BENCHMARK_TEMPLATE(BM_trsv_TransLower_tmpc_Static, double, 1, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_trsv_TransLower_tmpc_Static, double, 4, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_trsv_TransLower_tmpc_Static, double, 35, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_trsv_TransLower_tmpc_Static, double, 60, blaze::columnMajor);
}
