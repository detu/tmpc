#include <tmpc/math/Llh.hpp>

#include <blaze/Math.h>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Real, size_t M, size_t N>
    static void BM_SyrkPotrf_blaze_Static(::benchmark::State& state)
    {
        blaze::StaticMatrix<Real, M, N, blaze::columnMajor> A;
        blaze::StaticMatrix<Real, N, N, blaze::rowMajor> B;
        blaze::SymmetricMatrix<blaze::StaticMatrix<Real, N, N, blaze::rowMajor>> B0;

        randomize(A);
        makePositiveDefinite(B0);
        
        for (auto _ : state)
        {
            B = declsym(B0) + declsym(trans(A) * A);        
            tmpc::llh(B);
        }
    }

    
    BENCHMARK_TEMPLATE(BM_SyrkPotrf_blaze_Static, double, 4, 5);
    BENCHMARK_TEMPLATE(BM_SyrkPotrf_blaze_Static, double, 30, 35);
}
