#include <casadi_math_bench.h>

#include <tmpc/casadi_interface/GeneratedFunction.hpp>

#include <benchmark/benchmark.h>

#include <vector>


#define ADD_BM_MTIMES(mnpT) BENCHMARK_CAPTURE(BM_MTimes, mnpT, MTimes_ ## mnpT ## _functions())


namespace tmpc :: benchmark
{
    using casadi_interface::GeneratedFunction;


    static void BM_MTimes(::benchmark::State& state, casadi_functions const * fun)
    {
        GeneratedFunction f {fun};

        std::vector<casadi_real> const A(f.n_row_in(0) * f.n_col_in(0));
        std::vector<casadi_real> const B(f.n_row_in(1) * f.n_col_in(1));
        std::vector<casadi_real> C(f.n_row_out(0) * f.n_col_out(0));
        
        for (auto _ : state)
            f({A.data(), B.data()}, {C.data()});
    }


    ADD_BM_MTIMES(2x2x2_MX);
    ADD_BM_MTIMES(3x3x3_MX);
    ADD_BM_MTIMES(5x5x5_MX);
    ADD_BM_MTIMES(10x10x10_MX);
    ADD_BM_MTIMES(20x20x20_MX);
    ADD_BM_MTIMES(30x30x30_MX);
    ADD_BM_MTIMES(2x2x2_SX);
    ADD_BM_MTIMES(3x3x3_SX);
    ADD_BM_MTIMES(5x5x5_SX);
    ADD_BM_MTIMES(10x10x10_SX);
    ADD_BM_MTIMES(20x20x20_SX);
    ADD_BM_MTIMES(30x30x30_SX);
}
