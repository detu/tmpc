#include <blaze/Math.h>
#include <Eigen/Dense>

#include <tmpc/math/Trsv.hpp>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Real>
    static void BM_trsv_Blaze_Dynamic(::benchmark::State& state)
    {
        size_t const M = state.range(0);
        blaze::LowerMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>> A(M, M);        
        blaze::DynamicVector<Real, blaze::columnVector> B(M);
        blaze::DynamicVector<Real, blaze::columnVector> C(M);

        randomize(A);
        randomize(B);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(C = inv(A) * B);
    }


    template <typename Real, size_t M>
    static void BM_trsv_Blaze_Static(::benchmark::State& state)
    {
        blaze::LowerMatrix<blaze::StaticMatrix<Real, M, M, blaze::columnMajor>> A;        
        blaze::StaticVector<Real, M, blaze::columnVector> B;
        blaze::StaticVector<Real, M, blaze::columnVector> C;

        randomize(A);
        randomize(B);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(C = inv(A) * B);
    }


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


#if BLAZE_BLAS_MODE
    template <typename Real, bool SO>
    static void BM_trsv_Blaze_blas(::benchmark::State& state)
    {
        size_t const M = state.range(0);
        blaze::DynamicMatrix<Real, SO> A(M, M);        
        blaze::DynamicVector<Real, SO> B(M);
        blaze::DynamicVector<Real, SO> C(M);

        randomize(A);
        randomize(B);
        
        for (auto _ : state)
        {
            C = B;
            trsv(A, C, 'L', 'N', 'N');
        }
    }
#endif


    static void trsvBenchArguments(::benchmark::internal::Benchmark* b) 
    {
        b->Arg(1)->Arg(4)->Arg(35)->Arg(60);
    }


    BENCHMARK_TEMPLATE(BM_trsv_Blaze_Dynamic, double)->Apply(trsvBenchArguments);
    
    BENCHMARK_TEMPLATE(BM_trsv_Blaze_Static, double, 1);
    BENCHMARK_TEMPLATE(BM_trsv_Blaze_Static, double, 4);
    BENCHMARK_TEMPLATE(BM_trsv_Blaze_Static, double, 35);
    BENCHMARK_TEMPLATE(BM_trsv_Blaze_Static, double, 60);

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

#if BLAZE_BLAS_MODE
    BENCHMARK_TEMPLATE(BM_trsv_Blaze_blas, double, blaze::rowMajor)->Apply(trsvBenchArguments);
    BENCHMARK_TEMPLATE(BM_trsv_Blaze_blas, double, blaze::columnMajor)->Apply(trsvBenchArguments);
#endif
}
