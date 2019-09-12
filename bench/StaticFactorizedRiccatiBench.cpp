#include <tmpc/qp/MpipmWorkspace.hpp>
#include <tmpc/qp/StaticFactorizedRiccati.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>

#include "RiccatiBench.hpp"


namespace tmpc :: benchmark
{
    template <size_t NX, size_t NU>
    void BM_StaticFactorizedRiccati(::benchmark::State& state)
    {
        size_t const N = state.range(0);

        OcpGraph const g = ocpGraphLinear(N + 1);
        auto const sz = ocpSizeNominalMpc(N, NX, NU, 0, 0, 0, false);
        MpipmWorkspace<double> ws(g, sz);
        StaticFactorizedRiccati<double, NX, NU> riccati(g);

        randomizeQp(ws);

        // Disable openblas multithreading
        // openblas_set_num_threads(1);

        for (auto _ : state)
            riccati(ws, ws);
    }


    BENCHMARK_TEMPLATE(BM_StaticFactorizedRiccati, 4, 1)->Arg(100);
    BENCHMARK_TEMPLATE(BM_StaticFactorizedRiccati, 30, 5)->Arg(100);

    //BENCHMARK(BM_FactorizedRiccati)->Args({1, 2, 1});
}
