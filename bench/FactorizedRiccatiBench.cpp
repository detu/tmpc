#include <tmpc/qp/MpipmWorkspace.hpp>
#include <tmpc/qp/FactorizedRiccati.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>

#include "RiccatiBench.hpp"


namespace tmpc :: benchmark
{
    void BM_FactorizedRiccati(::benchmark::State& state)
    {
        size_t const N = state.range(0), nx = state.range(1), nu = state.range(2);

        OcpGraph const g = ocpGraphLinear(N + 1);
        auto const sz = ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true);
        MpipmWorkspace<double> ws(g, sz);
        FactorizedRiccati<double> riccati(g, sz);

        randomizeQp(ws);

        // Disable openblas multithreading
        // openblas_set_num_threads(1);

        for (auto _ : state)
            riccati(ws, ws);
    }


    BENCHMARK(BM_FactorizedRiccati)->Apply(riccatiBenchArguments);

    //BENCHMARK(BM_FactorizedRiccati)->Args({1, 2, 1});
}
