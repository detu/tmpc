#include <tmpc/qp/MpipmWorkspace.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>

#include "RiccatiBench.hpp"


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


    BENCHMARK(BM_MpipmRiccati)->Apply(riccatiBenchArguments);

    //BENCHMARK(BM_MpipmRiccati)->Args({100, 10, 5});
}
