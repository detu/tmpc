#include <blaze/Math.h>
#include <Eigen/Dense>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Real>
    static void BM_trmv_Blaze_Dynamic(::benchmark::State& state)
    {
        size_t const M = state.range(0);
        blaze::LowerMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>> A(M, M);        
        blaze::DynamicVector<Real, blaze::columnVector> B(M);
        blaze::DynamicVector<Real, blaze::columnVector> C(M);

        randomize(A);
        randomize(B);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(C = trans(A) * B);
    }


    template <typename Real, size_t M>
    static void BM_trmv_Blaze_Static(::benchmark::State& state)
    {
        blaze::LowerMatrix<blaze::StaticMatrix<Real, M, M, blaze::columnMajor>> A;        
        blaze::StaticVector<Real, M, blaze::columnVector> B;
        blaze::StaticVector<Real, M, blaze::columnVector> C;

        randomize(A);
        randomize(B);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(C = trans(A) * B);
    }


#if BLAZE_BLAS_MODE
    template <typename Real>
    static void BM_trmv_Blaze_blas(::benchmark::State& state)
    {
        size_t const M = state.range(0);
        blaze::DynamicMatrix<Real, blaze::columnMajor> A(M, M);        
        blaze::DynamicVector<Real, blaze::columnVector> B(M);
        blaze::DynamicVector<Real, blaze::columnVector> C(M);

        randomize(A);
        randomize(B);
        
        for (auto _ : state)
        {
            C = B;
            trmv(C, trans(std::as_const(A)), CblasUpper);
        }
    }
#endif


    static void trmvBenchArguments(::benchmark::internal::Benchmark* b) 
    {
        b->Arg(1)->Arg(4)->Arg(35);
    }


    BENCHMARK_TEMPLATE(BM_trmv_Blaze_Dynamic, double)->Apply(trmvBenchArguments);
    BENCHMARK_TEMPLATE(BM_trmv_Blaze_Static, double, 1);
    BENCHMARK_TEMPLATE(BM_trmv_Blaze_Static, double, 4);
    BENCHMARK_TEMPLATE(BM_trmv_Blaze_Static, double, 35);

#if BLAZE_BLAS_MODE
    BENCHMARK_TEMPLATE(BM_trmv_Blaze_blas, double)->Apply(trmvBenchArguments);
#endif
}
