#include <tmpc/qp/StaticOcpQp.hpp>
#include <tmpc/ocp/StaticOcpSolution.hpp>
#include <tmpc/qp/StaticClassicalRiccati.hpp>
#include <tmpc/qp/Randomize.hpp>

#include <bench/Complexity.hpp>


namespace tmpc :: benchmark
{
    template <size_t NX, size_t NU>
    void BM_StaticClassicalRiccati(State& state)
    {
        size_t const N = state.range(0);

        OcpTree const g(N + 1);
        StaticOcpQp<double, NX, NU> qp(g);
        StaticOcpSolution<double, NX, NU> sol(g);
        StaticClassicalRiccati<double, NX, NU> riccati(g);

        randomize(qp);

        for (auto _ : state)
            riccati(qp, sol);

        setCounters(state.counters, complexityClassicalRiccati(NX, NU, N));
        state.counters["nx"] = NX;
        state.counters["nu"] = NU;
        state.counters["n"] = N;
    }


#define BOOST_PP_LOCAL_LIMITS (1, BENCHMARK_RICCATI_MAX)
#define BOOST_PP_LOCAL_MACRO(n) \
    BENCHMARK_TEMPLATE(BM_StaticClassicalRiccati, n, n)->Arg(100);
#include BOOST_PP_LOCAL_ITERATE()
}
