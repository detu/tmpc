#include <tmpc/ocp/OcpGraph.hpp>

#include <benchmark/benchmark.h>

namespace tmpc :: benchmark
{
    class BfsBenchmarkVisitor
    :   public graph::default_bfs_visitor 
    {
    public:
        BfsBenchmarkVisitor(size_t& vertex_count, size_t& edge_count)
        :   vertexCount_(vertex_count)
        ,   edgeCount_(edge_count)
        {
            vertexCount_ = 0;
            edgeCount_ = 0;
        }

    
        template <typename Vertex, typename Graph>
        void discover_vertex(Vertex u, const Graph& g)
        {
            ++vertexCount_;
        }


        template <typename Edge, typename Graph>
        void tree_edge(Edge e, const Graph& g)
        {
            ++edgeCount_;
        }

    
    private:
        size_t vertexCount_;
        size_t edgeCount_;
    };


    class DfsBenchmarkVisitor
    :   public graph::default_dfs_visitor 
    {
    public:
        DfsBenchmarkVisitor(size_t& vertex_count, size_t& edge_count)
        :   vertexCount_(vertex_count)
        ,   edgeCount_(edge_count)
        {
            vertexCount_ = 0;
            edgeCount_ = 0;
        }

    
        template <typename Vertex, typename Graph>
        void discover_vertex(Vertex u, const Graph& g)
        {
            ++vertexCount_;
        }


        template <typename Edge, typename Graph>
        void tree_edge(Edge e, const Graph& g)
        {
            ++edgeCount_;
        }

    
    private:
        size_t vertexCount_;
        size_t edgeCount_;
    };


    void BM_OcpGraph_BreadthFirstSearch(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpGraph const g = ocpGraphLinear(N + 1);
        size_t nv = 0, ne = 0;

        for (auto _ : state)
            breadth_first_search(g, vertex(0, g), visitor(BfsBenchmarkVisitor(nv, ne)));
    }


    void BM_OcpGraph_DepthFirstSearch(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpGraph const g = ocpGraphLinear(N + 1);
        size_t nv = 0, ne = 0;

        for (auto _ : state)
            depth_first_search(g, visitor(DfsBenchmarkVisitor(nv, ne)));
    }


    void BM_OcpGraph_EnumVertices(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpGraph const g = ocpGraphLinear(N + 1);
        size_t nv = 0;

        for (auto _ : state)
            for(auto v : graph::vertices(g))
                ++nv;
    }


    void BM_OcpGraph_EnumEdges(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpGraph const g = ocpGraphLinear(N + 1);
        size_t ne = 0;

        for (auto _ : state)
            for(auto e : graph::edges(g))
                ++ne;
    }


    void BM_OcpGraph_Source(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpGraph const g = ocpGraphLinear(N + 1);
        OcpVertexDescriptor u;

        for (auto _ : state)
            for(auto e : graph::edges(g))
                u = source(e, g);
    }


    void BM_OcpGraph_Target(::benchmark::State& state)
    {
        size_t const N = 100;

        OcpGraph const g = ocpGraphLinear(N + 1);
        OcpVertexDescriptor v;

        for (auto _ : state)
            for(auto e : graph::edges(g))
                v = target(e, g);
    }


    BENCHMARK(BM_OcpGraph_BreadthFirstSearch);
    BENCHMARK(BM_OcpGraph_DepthFirstSearch);
    BENCHMARK(BM_OcpGraph_EnumVertices);
    BENCHMARK(BM_OcpGraph_EnumEdges);
    BENCHMARK(BM_OcpGraph_Source);
    BENCHMARK(BM_OcpGraph_Target);
}
