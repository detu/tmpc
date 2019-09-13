#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/graph/Graph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/graph/DepthFirstSearch.hpp>
#include <tmpc/graph/ImpactRecorder.hpp>
#include <tmpc/core/PropertyMap.hpp>

#include <boost/iterator/iterator_facade.hpp>

#include <vector>
#include <optional>


namespace tmpc
{
    class OcpGraph
    :   public graph::compressed_sparse_row_graph<graph::directedS>
    {
    public:
        OcpGraph() = default;


        template <typename InputIterator>
        OcpGraph(InputIterator first, InputIterator last, vertices_size_type numverts)
        :   Base(boost::edges_are_sorted, first, last, numverts)
        {
            auto& base = static_cast<Base const&>(*this);
            
            parentEdge_.resize(num_vertices(base));
            for (auto e : graph::edges(base))
                parentEdge_[target(e, base)] = e;
        }


        std::optional<edge_descriptor> parentEdge(vertex_descriptor v) const
        {
            std::optional<edge_descriptor> e;

            if (v != 0)
                e = parentEdge_[v];

            return e;
        } 


    private:
        using Base = boost::compressed_sparse_row_graph<boost::directedS>;

        // Maps vertex id to its parent edge.
        std::vector<edge_descriptor> parentEdge_;
    };

    // using OcpGraph = boost::adjacency_list<
    //     boost::vecS, boost::vecS, boost::bidirectionalS,
    //     boost::no_property, boost::property<boost::edge_index_t, size_t>
    // >;


    using OcpVertexDescriptor = boost::graph_traits<OcpGraph>::vertex_descriptor;
    using OcpEdgeDescriptor = boost::graph_traits<OcpGraph>::edge_descriptor;


    namespace detail
    {
        template <typename OutDegreeIterator>
        class OutDegreeListToEdgesIterator
        :   public boost::iterator_facade<
            OutDegreeListToEdgesIterator<OutDegreeIterator>,
            std::pair<size_t, size_t> const,  // Value
            std::forward_iterator_tag  // CategoryOrTraversal
            //std::pair<size_t, size_t>   // Reference
        >
        {
        public:
            OutDegreeListToEdgesIterator()
            :   uv_(0, 0)
            ,   k_(0)
            {
            }


            OutDegreeListToEdgesIterator(OutDegreeListToEdgesIterator const&) = default;


            explicit OutDegreeListToEdgesIterator(OutDegreeIterator out_degree)
            :   od_(out_degree)
            ,   uv_(0, 1)
            ,   k_(*od_)
            {
            }


        private:
            friend class boost::iterator_core_access;


            std::pair<size_t, size_t> const& dereference() const
            {
                if (k_ == 0)
                    throw std::logic_error("Dereferencing an end-of-sequence OutDegreeListToEdgesIterator");

                return uv_;
            }


            void increment() 
            { 
                if (k_ == 0)
                    throw std::logic_error("Incrementing an end-of-sequence OutDegreeListToEdgesIterator");

                --k_;

                while (k_ == 0 && uv_.first < uv_.second)
                {
                    k_ = *++od_;
                    ++uv_.first;
                }

                if (k_ > 0)
                    ++uv_.second;                
            }


            bool equal(OutDegreeListToEdgesIterator const& other) const
            {
                return (od_ == other.od_ && uv_ == other.uv_) || (k_ == 0 && other.k_ == 0);
            }


            OutDegreeIterator od_;
            std::pair<size_t, size_t> uv_;
            size_t k_;
        };
    }


    /// Create a tree-structured OcpGraph from a list defined by the out_gedree iterator,
    /// which is the number out-edges of each node, in breadth-first order.
    ///
    template <typename InputIterator>
    inline OcpGraph ocpGraphFromOutDegree(InputIterator out_degree)
    {
        detail::OutDegreeListToEdgesIterator<InputIterator> first(out_degree), last;
        size_t const n_edges = std::distance(first, last);
        
        return OcpGraph(first, last, n_edges + 1);
    }


    /// @brief Create a linear OcpGraph with n_nodes nodes.
    ///
    OcpGraph ocpGraphLinear(size_t n_nodes);


    /// @brief Create a robust MPC OcpGraph.
    ///
    /// @param depth depth of the graph. It is equal to number of control intervals + 1.
    /// depth = 0 means no nodes; depth = 1 means one layer of nodes (consisting of the root node only);
    /// depth = 2 means the root node and 1 layer of nodes after it (1 time step) and so on.
    /// @param branching how many children each parent node has within the robust horizon
    /// @param robust_horizon time step after which branching stops. robust_horizon = 0 means
    /// no branching, robust_horizon = 1 means branching only at the root node.
    OcpGraph ocpGraphRobustMpc(size_t depth, size_t branching, size_t robust_horizon);


    /// Vertex index property map.
    inline auto vertexIndex(OcpGraph const& g)
    {
        return get(graph::vertex_index, g);
    }


    /// @brief Input degree of a node.
    ///
    /// Since OcpGraph is always a tree, in_degree() is 0 for the root and 1 for all other nodes.
    inline OcpGraph::edges_size_type in_degree(OcpVertexDescriptor v, OcpGraph const& g)
    {
        return get(graph::vertex_index, g, v) == 0 ? 0 : 1;
    }
	

    /// @brief Parent of a node.
    ///
    /// Since OcpGraph is always a tree, parent() is empty for the root 
    /// and equals to source(in_edges(v, g), g) for non-root nodes.
    inline std::optional<OcpVertexDescriptor> parent(OcpVertexDescriptor v, OcpGraph const& g)
    {
        std::optional<OcpVertexDescriptor> p;

        if (auto const e = g.parentEdge(v))
            p = source(*e, g);

        return p;
    }


    /// @brief Siblings of a vertex.
    ///
    /// Siblings of v is the set of vertices u such as parent(u, g) == parent(v, g).
    inline auto siblings(OcpVertexDescriptor v, OcpGraph const& g)
    {
        auto const p = parent(v, g);

        decltype(graph::adjacent_vertices(*p, g)) sb;

        if (p)
            sb = graph::adjacent_vertices(*p, g);

        return sb;
    }


    /// @brief Root of the graph.
    inline OcpVertexDescriptor root(OcpGraph const& g)
    {
        return vertex(0, g);
    }


    // template <typename Graph, typename ImpactMap, typename ColorMap>
    // inline void recordImpact(Graph const& g, ImpactMap impact, ColorMap color)
    // {
    //     depth_first_visit(g, root(g), dfs_visitor(ImpactRecorder(impact)), color);
    // }


    // template <typename Graph, typename ImpactMap>
    // inline void recordImpact(Graph const& g, ImpactMap impact)
    // {
    //     std::vector<boost::default_color_type> color(num_vertices(g));

    //     depth_first_visit(g, root(g), dfs_visitor(ImpactRecorder(impact)), 
    //         make_iterator_property_map(color.begin(), get(graph::vertex_index, g)));
    // }


    /// @brief Traverse the OcpGraph in breadth-first order (from root to leaves)
    template <typename Visitor>
    inline void breadthFirstVisit(OcpGraph const& g, Visitor visitor)
    {
        for (auto v : graph::vertices(g))
        {
            if (auto const e = g.parentEdge(v))
                visitor.edge(*e, g);

            visitor.vertex(v, g);
        }
    }


    /// @brief Traverse the OcpGraph in reverse breadth-first order (from leaves to root)
    template <typename Visitor>
    inline void reverseBreadthFirstVisit(OcpGraph const& g, Visitor visitor)
    {
        // Beware of Incorrect behavior of reversed_range with -O2 g++ option:
		// https://github.com/boostorg/range/issues/82
        // for (auto v : reverse(graph::vertices(g)))

        for (OcpVertexDescriptor v = num_vertices(g); v-- > 0; )
        {
            visitor.vertex(v, g);

            if (auto const e = g.parentEdge(v))
                visitor.edge(*e, g);
        }
    }
}
