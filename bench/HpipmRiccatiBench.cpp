#include <tmpc/qp/HpipmWorkspace.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    void BM_HpipmRiccati(::benchmark::State& state)
    {
        size_t const N = state.range(0), nx = state.range(1), nu = state.range(2);

        OcpGraph const g = ocpGraphLinear(N + 1);
        HpipmWorkspace<BlazeKernel<double>> ws(g, ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true));

        randomizeQp(ws);
        ws.convertProblem();

        for (auto _ : state)
            ws.solveUnconstrainedInternal();
    }

    // static void BenchmarkArguments(benchmark::internal::Benchmark* b) 
    // {
    //     for (int nx : {2, 5, 10})
    //         for (int nu : {1, 2})
    //             b->Args({100, nx, nu});
    // }


    
    // BENCHMARK(BM_HpipmRiccati)->Apply(BenchmarkArguments);


    BENCHMARK(BM_HpipmRiccati)->Args({100, 10, 5});
}
