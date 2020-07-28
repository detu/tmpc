#include <tmpc/qp/StaticOcpQp.hpp>
#include <tmpc/ocp/StaticOcpSolution.hpp>
#include <tmpc/qp/StaticFactorizedRiccati.hpp>
#include <tmpc/qp/Randomize.hpp>

#include <bench/Benchmark.hpp>
#include <bench/Complexity.hpp>

#include <blasfeo/Blasfeo.hpp>


// #define BENCHMARK_STATIC_FACTORIZED_RICCATI(NX) \
//     BENCHMARK_TEMPLATE(BM_StaticFactorizedRiccati, NX, 1)->Arg(100); \
//     BENCHMARK_TEMPLATE(BM_StaticFactorizedRiccati, NX, NX)->Arg(100);

#define BENCHMARK_STATIC_FACTORIZED_RICCATI(NX) \
    BENCHMARK_TEMPLATE(BM_StaticFactorizedRiccati, NX, NX)->Arg(100);


namespace tmpc :: benchmark
{
    template <size_t NX, size_t NU>
    void BM_StaticFactorizedRiccati(State& state)
    {
        size_t const N = state.range(0);

        OcpTree const g(N + 1);
        StaticOcpQp<double, NX, NU> qp(g);
        StaticOcpSolution<double, NX, NU> sol(g);
        StaticFactorizedRiccati<double, NX, NU> riccati(g);

        randomize(qp);

        for (auto _ : state)
            riccati(qp, sol);

        setCounters(state.counters, complexityFactorizedRiccati(NX, NU, N));
        state.counters["nx"] = NX;
        state.counters["nu"] = NU;
        state.counters["n"] = N;
    }


#define BOOST_PP_LOCAL_LIMITS (1, 20)
#define BOOST_PP_LOCAL_MACRO(n) BENCHMARK_STATIC_FACTORIZED_RICCATI(n);
#include BOOST_PP_LOCAL_ITERATE()
}
