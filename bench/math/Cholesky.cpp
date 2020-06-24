#include <tmpc/math/Llh.hpp>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Real>
    static void BM_Cholesky_tmpc_llh_Dynamic(::benchmark::State& state)
    {
        size_t const N = state.range(0);
        blaze::DynamicMatrix<Real, blaze::columnMajor> A(N, N);        
        makePositiveDefinite(A);

        blaze::DynamicMatrix<Real, blaze::columnMajor> B(N, N);        
        for (auto _ : state)
            tmpc::llh(A, B);
    }


    template <typename Real, int N, bool SO1, bool SO2>
    static void BM_Cholesky_tmpc_llh_Static(::benchmark::State& state)
    {
        // Make a positive definite matrix A
        blaze::StaticMatrix<Real, N, N, SO1> A;
        makePositiveDefinite(A);

        blaze::LowerMatrix<blaze::StaticMatrix<Real, N, N, SO2>> L;

        for (auto _ : state)
            tmpc::llh(A, L);
    }


    template <typename Real, int N, bool SO1, bool SO2>
    static void BM_Cholesky_tmpc_llh_Static_Symmetric(::benchmark::State& state)
    {
        // Make a positive definite matrix A
        blaze::SymmetricMatrix<blaze::StaticMatrix<Real, N, N, SO1>> A;
        makePositiveDefinite(A);

        blaze::LowerMatrix<blaze::StaticMatrix<Real, N, N, SO2>> L;

        for (auto _ : state)
            tmpc::llh(A, L);
    }


    template <typename Real, int N, bool SO>
    static void BM_Cholesky_tmpc_llh_Static_InPlace(::benchmark::State& state)
    {
        // Make a positive definite matrix A
        blaze::StaticMatrix<Real, N, N, SO> A, A0;
        makePositiveDefinite(A0);

        for (auto _ : state)
        {
            state.PauseTiming();
            A = A0;
            state.ResumeTiming();

            tmpc::llh(A);
        }
    }


    static void choleskyBenchArguments(::benchmark::internal::Benchmark* b) 
    {
        b->Arg(1)->Arg(2)->Arg(5)->Arg(10)->Arg(35);
    }

    
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Dynamic, double)->Apply(choleskyBenchArguments);
    
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 1, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 2, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 5, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 10, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 35, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 1, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 2, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 5, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 10, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 35, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 1, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 2, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 5, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 10, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static, double, 35, blaze::rowMajor, blaze::rowMajor);

    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 1, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 2, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 5, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 10, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 35, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 1, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 2, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 5, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 10, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_InPlace, double, 35, blaze::rowMajor);

    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 1, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 2, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 5, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 10, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 35, blaze::columnMajor, blaze::columnMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 1, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 2, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 5, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 10, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 35, blaze::columnMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 1, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 2, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 5, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 10, blaze::rowMajor, blaze::rowMajor);
    BENCHMARK_TEMPLATE(BM_Cholesky_tmpc_llh_Static_Symmetric, double, 35, blaze::rowMajor, blaze::rowMajor);
}
