#include <tmpc/math/SyrkPotrf.hpp>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Real, size_t M, size_t N>
    static void BM_SyrkPotrf_tmpc_Static(::benchmark::State& state)
    {
        blaze::StaticMatrix<Real, M, N, blaze::columnMajor> A;
        blaze::SymmetricMatrix<blaze::StaticMatrix<Real, N, N, blaze::columnMajor>> C;
        blaze::LowerMatrix<blaze::StaticMatrix<Real, N, N, blaze::columnMajor>> D;

        randomize(A);
        makePositiveDefinite(C);
        
        for (auto _ : state)
            tmpc::syrkPotrf(A, declsym(C), D);
    }


    BENCHMARK_TEMPLATE(BM_SyrkPotrf_tmpc_Static, double, 4, 5);
    BENCHMARK_TEMPLATE(BM_SyrkPotrf_tmpc_Static, double, 30, 60);
}
