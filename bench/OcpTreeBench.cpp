#include <tmpc/ocp/OcpTree.hpp>

#include <benchmark/benchmark.h>

namespace tmpc :: benchmark
{
    void BM_OcpTree_EnumVertices(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpTree const g(N + 1);
        size_t nv = 0;

        for (auto _ : state)
            for(auto v : vertices(g))
                ++nv;
    }


    void BM_OcpTree_EnumEdges(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpTree const g(N + 1);
        size_t ne = 0;

        for (auto _ : state)
            for(auto e : edges(g))
                ++ne;
    }


    void BM_OcpTree_Source(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpTree const g(N + 1);
        OcpVertex u;

        for (auto _ : state)
            for(auto e : edges(g))
                u = source(e, g);
    }


    void BM_OcpTree_Target(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpTree const g(N + 1);
        OcpVertex v;

        for (auto _ : state)
            for(auto e : edges(g))
                v = target(e, g);
    }


    BENCHMARK(BM_OcpTree_EnumVertices);
    BENCHMARK(BM_OcpTree_EnumEdges);
    BENCHMARK(BM_OcpTree_Source);
    BENCHMARK(BM_OcpTree_Target);
}
