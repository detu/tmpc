#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/qp/OcpQp.hpp>

#include <benchmark/benchmark.h>


namespace tmpc :: benchmark
{
    template <typename Workspace>
    void BM_Riccati(::benchmark::State& state)
    {
        size_t const N = state.range(0), nx = state.range(1), nu = state.range(2);

        OcpGraph const g = ocpGraphLinear(N + 1);
        Workspace ws(g, ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true));

        randomizeQp(ws);

        //blasfeo_pack_dmat(int m, int n, double *A, int lda, struct blasfeo_dmat *sB, int bi, int bj);
        
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


    //BENCHMARK(BM_MpipmRiccati)->Args({100, 10, 5});
}
