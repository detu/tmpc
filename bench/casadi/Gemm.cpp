#include <generated_gemm.h>

#include <tmpc/casadi/GeneratedFunction.hpp>

#include <benchmark/benchmark.h>

#include <vector>


#define ADD_BM_GEMM(mnpT) BENCHMARK_CAPTURE(BM_gemm, mnpT, generated_gemm_ ## mnpT ## _functions())


namespace tmpc :: benchmark
{
    static void BM_gemm(::benchmark::State& state, casadi_functions const * fun)
    {
        casadi::GeneratedFunction f {fun};

        std::vector<casadi_real> const A(f.n_row_in(0) * f.n_col_in(0));
        std::vector<casadi_real> const B(f.n_row_in(1) * f.n_col_in(1));
        std::vector<casadi_real> const C(f.n_row_in(2) * f.n_col_in(2));
        std::vector<casadi_real> D(f.n_row_out(0) * f.n_col_out(0));
        
        for (auto _ : state)
            f({A.data(), B.data(), C.data()}, {D.data()});
    }


    ADD_BM_GEMM(2x2x2_MX);
    ADD_BM_GEMM(3x3x3_MX);
    ADD_BM_GEMM(5x5x5_MX);
    ADD_BM_GEMM(10x10x10_MX);
    ADD_BM_GEMM(20x20x20_MX);
    ADD_BM_GEMM(30x30x30_MX);
    ADD_BM_GEMM(2x2x2_SX);
    ADD_BM_GEMM(3x3x3_SX);
    ADD_BM_GEMM(5x5x5_SX);
    ADD_BM_GEMM(10x10x10_SX);
    ADD_BM_GEMM(20x20x20_SX);
    ADD_BM_GEMM(30x30x30_SX);
}
