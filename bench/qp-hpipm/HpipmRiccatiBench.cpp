#include <tmpc/hpipm/NominalSolver.hpp>
#include <tmpc/qp/Randomize.hpp>
#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/ocp/DynamicOcpSolution.hpp>

#include <bench/qp/DynamicRiccatiBench.hpp>
#include <bench/Benchmark.hpp>
#include <bench/Complexity.hpp>


namespace tmpc :: benchmark
{
    void BM_HpipmRiccati(State& state)
    {
        size_t const N = state.range(0), nx = state.range(1), nu = state.range(2);

        auto const size = ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true);
        DynamicOcpQp<double> qp {size};
        DynamicOcpSolution<double> sol {size};
        NominalSolver<double> ws {size};

        randomize(qp);
        ws.convertProblem(qp);

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
