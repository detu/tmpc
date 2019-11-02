#include <tmpc/blasfeo/Blasfeo.hpp>

#include <blaze/Math.h>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    void BM_potrf(::benchmark::State& state)
    {
        size_t const m = state.range(0);

        // Init Blaze matrices
        //
        blaze::DynamicMatrix<double, blaze::columnMajor> blaze_A(m, m), blaze_L(m, m);
        makePositiveDefinite(blaze_A);
        
        // Init BLASFEO matrices
        //
        blasfeo::DynamicMatrix<double> blasfeo_A(blaze_A), blasfeo_L(m, m);
        
        // Do potrf with BLASFEO
        for (auto _ : state)
            potrf(m, blasfeo_A, 0, 0, blasfeo_L, 0, 0);
    }


    BENCHMARK(BM_potrf)->DenseRange(1, 50);
}
