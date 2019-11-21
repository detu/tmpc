#include <tmpc/qp/HpipmWorkspace.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <bench/RiccatiBench.hpp>
#include <bench/Benchmark.hpp>
#include <bench/Complexity.hpp>


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

        setCounters(state.counters, complexityFactorizedRiccati(nx, nu, N));
        state.counters["nx"] = nx;
        state.counters["nu"] = nu;
        state.counters["n"] = N;
    }


    BENCHMARK(BM_HpipmRiccati)->Apply(riccatiBenchArguments);

    //BENCHMARK(BM_HpipmRiccati)->Args({100, 10, 5});
}
