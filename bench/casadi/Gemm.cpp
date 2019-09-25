#include <casadi/casadi.hpp>    // must be included before <casadi/mem.h>, indirectly included from <generated_gemm.h>!
#include <generated_gemm.h>

#include <tmpc/casadi/GeneratedFunction.hpp>

#include <benchmark/benchmark.h>

#include <vector>


#define ADD_BM_GEMM_GENERATED(mnpT) BENCHMARK_CAPTURE(BM_gemm_generated, mnpT, generated_gemm_ ## mnpT ## _functions())
#define ADD_BM_GEMM_GENERATED1(m, T) BENCHMARK_CAPTURE(BM_gemm_generated, m/T, generated_gemm_ ## m ## x ## m ## x ## m ## _ ## T ## _functions())


namespace tmpc :: benchmark
{
    using namespace ::benchmark;


    static void BM_gemm_generated(::benchmark::State& state, casadi_functions const * fun)
    {
        casadi::GeneratedFunction f {fun};

        std::vector<casadi_real> const A(f.n_row_in(0) * f.n_col_in(0));
        std::vector<casadi_real> const B(f.n_row_in(1) * f.n_col_in(1));
        std::vector<casadi_real> const C(f.n_row_in(2) * f.n_col_in(2));
        std::vector<casadi_real> D(f.n_row_out(0) * f.n_col_out(0));
        
        for (auto _ : state)
            f({A.data(), B.data(), C.data()}, {D.data()});

        state.counters["flops"] = Counter(f.n_row_in(0) * f.n_col_in(0) * f.n_col_in(1), Counter::kIsIterationInvariantRate);
        state.counters["m"] = f.n_row_in(0);
    }


    template <typename X>
    static void BM_gemm_native(::benchmark::State& state)
    {
        size_t const m = state.range(0);

        X A = X::sym("a", m, m);
        X B = X::sym("b", m, m);
        X C = X::sym("c", m, m);

        ::casadi::Function f("f", std::vector<X> {A, B, C}, std::vector<X> {C + mtimes(transpose(A), B)});

        std::vector<double> const a(m * m);
        std::vector<double> const b(m * m);
        std::vector<double> const c(m * m);
        std::vector<double> d(m * m);

        std::vector<double const *> in {a.data(), b.data(), c.data()};
        std::vector<double *> out {d.data()};

        for (auto _ : state)
            f(in, out);

        state.counters["flops"] = Counter(m * m * m, Counter::kIsIterationInvariantRate);
        state.counters["m"] = m;
    }


    ADD_BM_GEMM_GENERATED1(1, SX);
    ADD_BM_GEMM_GENERATED1(2, SX);
    ADD_BM_GEMM_GENERATED1(3, SX);
    ADD_BM_GEMM_GENERATED1(4, SX);
    ADD_BM_GEMM_GENERATED1(5, SX);
    ADD_BM_GEMM_GENERATED1(6, SX);
    ADD_BM_GEMM_GENERATED1(7, SX);
    ADD_BM_GEMM_GENERATED1(8, SX);
    ADD_BM_GEMM_GENERATED1(9, SX);
    ADD_BM_GEMM_GENERATED1(10, SX);
    ADD_BM_GEMM_GENERATED1(11, SX);
    ADD_BM_GEMM_GENERATED1(12, SX);
    ADD_BM_GEMM_GENERATED1(13, SX);
    ADD_BM_GEMM_GENERATED1(14, SX);
    ADD_BM_GEMM_GENERATED1(15, SX);
    ADD_BM_GEMM_GENERATED1(16, SX);
    ADD_BM_GEMM_GENERATED1(17, SX);
    ADD_BM_GEMM_GENERATED1(18, SX);
    ADD_BM_GEMM_GENERATED1(19, SX);
    ADD_BM_GEMM_GENERATED1(20, SX);
    ADD_BM_GEMM_GENERATED1(21, SX);
    ADD_BM_GEMM_GENERATED1(22, SX);
    ADD_BM_GEMM_GENERATED1(23, SX);
    ADD_BM_GEMM_GENERATED1(24, SX);
    ADD_BM_GEMM_GENERATED1(25, SX);
    ADD_BM_GEMM_GENERATED1(26, SX);
    ADD_BM_GEMM_GENERATED1(27, SX);
    ADD_BM_GEMM_GENERATED1(28, SX);
    ADD_BM_GEMM_GENERATED1(29, SX);
    ADD_BM_GEMM_GENERATED1(30, SX);
    ADD_BM_GEMM_GENERATED1(31, SX);
    ADD_BM_GEMM_GENERATED1(32, SX);
    ADD_BM_GEMM_GENERATED1(33, SX);
    ADD_BM_GEMM_GENERATED1(34, SX);
    ADD_BM_GEMM_GENERATED1(35, SX);
    ADD_BM_GEMM_GENERATED1(36, SX);
    ADD_BM_GEMM_GENERATED1(37, SX);
    ADD_BM_GEMM_GENERATED1(38, SX);
    ADD_BM_GEMM_GENERATED1(39, SX);
    ADD_BM_GEMM_GENERATED1(40, SX);

    ADD_BM_GEMM_GENERATED1(1, MX);
    ADD_BM_GEMM_GENERATED1(2, MX);
    ADD_BM_GEMM_GENERATED1(3, MX);
    ADD_BM_GEMM_GENERATED1(4, MX);
    ADD_BM_GEMM_GENERATED1(5, MX);
    ADD_BM_GEMM_GENERATED1(6, MX);
    ADD_BM_GEMM_GENERATED1(7, MX);
    ADD_BM_GEMM_GENERATED1(8, MX);
    ADD_BM_GEMM_GENERATED1(9, MX);
    ADD_BM_GEMM_GENERATED1(10, MX);
    ADD_BM_GEMM_GENERATED1(11, MX);
    ADD_BM_GEMM_GENERATED1(12, MX);
    ADD_BM_GEMM_GENERATED1(13, MX);
    ADD_BM_GEMM_GENERATED1(14, MX);
    ADD_BM_GEMM_GENERATED1(15, MX);
    ADD_BM_GEMM_GENERATED1(16, MX);
    ADD_BM_GEMM_GENERATED1(17, MX);
    ADD_BM_GEMM_GENERATED1(18, MX);
    ADD_BM_GEMM_GENERATED1(19, MX);
    ADD_BM_GEMM_GENERATED1(20, MX);
    ADD_BM_GEMM_GENERATED1(21, MX);
    ADD_BM_GEMM_GENERATED1(22, MX);
    ADD_BM_GEMM_GENERATED1(23, MX);
    ADD_BM_GEMM_GENERATED1(24, MX);
    ADD_BM_GEMM_GENERATED1(25, MX);
    ADD_BM_GEMM_GENERATED1(26, MX);
    ADD_BM_GEMM_GENERATED1(27, MX);
    ADD_BM_GEMM_GENERATED1(28, MX);
    ADD_BM_GEMM_GENERATED1(29, MX);
    ADD_BM_GEMM_GENERATED1(30, MX);
    ADD_BM_GEMM_GENERATED1(31, MX);
    ADD_BM_GEMM_GENERATED1(32, MX);
    ADD_BM_GEMM_GENERATED1(33, MX);
    ADD_BM_GEMM_GENERATED1(34, MX);
    ADD_BM_GEMM_GENERATED1(35, MX);
    ADD_BM_GEMM_GENERATED1(36, MX);
    ADD_BM_GEMM_GENERATED1(37, MX);
    ADD_BM_GEMM_GENERATED1(38, MX);
    ADD_BM_GEMM_GENERATED1(39, MX);
    ADD_BM_GEMM_GENERATED1(40, MX);

    using ::casadi::SX;
    using ::casadi::MX;

    BENCHMARK_TEMPLATE(BM_gemm_native, SX)->DenseRange(1, 40);
    BENCHMARK_TEMPLATE(BM_gemm_native, MX)->DenseRange(1, 40);
}
