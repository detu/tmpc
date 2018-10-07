#include <tmpc/qp/MpipmWorkspace.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    void BM_MpipmRiccati(::benchmark::State& state)
    {
        size_t const N = state.range(0), nx = state.range(1), nu = state.range(2);

        OcpGraph const g = ocpGraphLinear(N + 1);
        MpipmWorkspace<double> ws(g, ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true));

        randomizeQp(ws);

        for (auto _ : state)
            ws.solveUnconstrained();
    }


    // static void BenchmarkArguments(benchmark::internal::Benchmark* b) 
    // {
    //     for (int nx : {2, 5, 10})
    //         for (int nu : {1, 2})
    //             b->Args({100, nx, nu});
    // }


    
    // BENCHMARK(BM_HpipmRiccati)->Apply(BenchmarkArguments);


    BENCHMARK(BM_MpipmRiccati)->Args({100, 10, 5});
}
