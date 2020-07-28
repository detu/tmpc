#include <tmpc/hpipm/NominalSolver.hpp>
#include <tmpc/qp/Randomize.hpp>
#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/ocp/DynamicOcpSolution.hpp>

#include <bench/qp/DynamicRiccatiBench.hpp>
#include <bench/Benchmark.hpp>
#include <bench/Complexity.hpp>


namespace tmpc :: benchmark
{
    static void BM_RiccatiBackend(State& state, int alg)
    {
        size_t const N = state.range(0), nx = state.range(1), nu = state.range(2);

        DynamicOcpSize const size(N, nx, nu, 0, 0, 0, true);
        DynamicOcpQp<double> qp {size};
        DynamicOcpSolution<double> sol {size};

        // Create options and set Riccati algiorithm to factorized.
        hpipm::OcpQpIpmArg<double> options {hpipm::OcpQpDim<double> {size}, ROBUST};
        options.set_ric_alg(alg);

        // Create solver
        hpipm::NominalSolver<double> ws {size, std::move(options)};

        randomize(qp);
        ws.convertProblem(qp);

        for (auto _ : state)
            ws.solveUnconstrainedInternal();

        switch (alg)
        {
            case 0:
            {
                setCounters(state.counters, complexityClassicalRiccati(nx, nu, N));
                break;
            }
            case 1:
            {
                setCounters(state.counters, complexityFactorizedRiccati(nx, nu, N));
                break;
            }
        }

        state.counters["nx"] = nx;
        state.counters["nu"] = nu;
        state.counters["n"] = N;
    }


    static void BM_ClassicalRiccati(State& state)
    {
        BM_RiccatiBackend(state, 0);
    }


    static void BM_FactorizedRiccati(State& state)
    {
        BM_RiccatiBackend(state, 1);
    }


    BENCHMARK(BM_ClassicalRiccati)->Apply(riccatiBenchArguments);
    BENCHMARK(BM_FactorizedRiccati)->Apply(riccatiBenchArguments);

    //BENCHMARK(BM_FactorizedRiccati)->Args({100, 10, 5});
}
