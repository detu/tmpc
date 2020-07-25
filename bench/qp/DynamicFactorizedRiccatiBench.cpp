#include <tmpc/qp/DynamicFactorizedRiccati.hpp>
#include <tmpc/qp/Randomize.hpp>
#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/ocp/DynamicOcpSolution.hpp>

#include <bench/qp/DynamicRiccatiBench.hpp>
#include <bench/Complexity.hpp>


namespace tmpc :: benchmark
{
    void BM_DynamicFactorizedRiccati(State& state)
    {
        using Real = double;
        size_t const N = state.range(0), nx = state.range(1), nu = state.range(2);

        auto const sz = ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true);
        DynamicOcpQp<Real> qp {sz};
        DynamicOcpSolution<Real> sol {sz};
        DynamicFactorizedRiccati<Real> riccati {sz};

        randomize(qp);

        for (auto _ : state)
            riccati(qp, sol);

        setCounters(state.counters, complexityFactorizedRiccati(nx, nu, N));
        state.counters["nx"] = nx;
        state.counters["nu"] = nu;
        state.counters["n"] = N;
    }

    
    BENCHMARK(BM_DynamicFactorizedRiccati)->Apply(riccatiBenchArguments);
}
