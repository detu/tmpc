#include <tmpc/blasfeo/Blasfeo.hpp>

#include <blaze/Math.h>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    using namespace ::benchmark;


    template <typename T>
    void BM_potrf(::benchmark::State& state)
    {
        size_t const m = state.range(0);

        // Init Blaze matrices
        //
        blaze::DynamicMatrix<T, blaze::columnMajor> blaze_A(m, m), blaze_L(m, m);
        makePositiveDefinite(blaze_A);
        
        // Init BLASFEO matrices
        //
        blasfeo::DynamicMatrix<T> blasfeo_A(blaze_A), blasfeo_L(m, m);
        
        // Do potrf with BLASFEO
        for (auto _ : state)
            potrf(m, blasfeo_A, 0, 0, blasfeo_L, 0, 0);

        state.counters["flops"] = Counter(m * m * m / 3., Counter::kIsIterationInvariantRate);
        state.counters["m"] = m;
    }


    BENCHMARK_TEMPLATE(BM_potrf, double)->DenseRange(1, 300);
}
