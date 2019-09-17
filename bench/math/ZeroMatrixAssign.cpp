#include <tmpc/EigenKernel.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <benchmark/benchmark.h>

#include <vector>


namespace tmpc :: benchmark
{
    template <typename Real, size_t M, size_t N>
    static void BM_DynamicMatrixZeroMatrixAssign(::benchmark::State& state)
    {
        blaze::DynamicMatrix<Real> A(M, N);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(A = blaze::ZeroMatrix<Real>(M, N));
    }


    template <typename Real, size_t M, size_t N>
    static void BM_DynamicMatrixZeroAssign(::benchmark::State& state)
    {
        blaze::DynamicMatrix<Real> A(M, N);
        
        for (auto _ : state)
            ::benchmark::DoNotOptimize(A = Real {0});
    }


    BENCHMARK_TEMPLATE(BM_DynamicMatrixZeroMatrixAssign, double, 4, 1);
    BENCHMARK_TEMPLATE(BM_DynamicMatrixZeroAssign, double, 4, 1);
}
