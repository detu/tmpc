#include <tmpc/qp/CondensingN3.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/EigenKernel.hpp>

#include <benchmark/benchmark.h>

#include <vector>


namespace tmpc :: benchmark
{
    static size_t constexpr NX = 13;
    static size_t constexpr NU = 8;
    static size_t constexpr NC = 0;
    static size_t constexpr NS = 0;

    template <typename Kernel>
    static auto makeTestQp(size_t N, size_t nx, size_t nu, size_t nc, size_t ns)
    {
        std::vector<OcpQp<Kernel>> qp;
        qp.reserve(N);

        for (size_t k = 0; k < N; ++k)
            qp.emplace_back(OcpSize {nx, nu, nc, ns}, nx);

        return qp;
    }

    template <typename CondensingAlgorithm>
    static void BM_Condensing(::benchmark::State& state)
    {
        using Kernel = typename CondensingAlgorithm::Kernel;
        
        auto qp = makeTestQp<Kernel>(state.range(0), NX, NU, NC, NS);
        CondensingAlgorithm condensing {sizeBegin(qp), sizeEnd(qp)};

        OcpQp<Kernel> condensed {condensing.condensedSize(), rows(qp.back().A())};
        for (auto _ : state)
            condensed = condensing(qp.begin(), qp.end());
    }


    BENCHMARK_TEMPLATE(BM_Condensing, CondensingN3<BlazeKernel<double>>)->Arg(5)->Arg(8)->Range(10, 100);
    BENCHMARK_TEMPLATE(BM_Condensing, CondensingN3<EigenKernel<double>>)->Arg(5)->Arg(8)->Range(10, 100);
}


BENCHMARK_MAIN();
