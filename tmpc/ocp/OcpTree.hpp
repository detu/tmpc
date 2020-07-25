#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/Exception.hpp>
#include <tmpc/property_map/PropertyMap.hpp>

#include <vector>
#include <optional>
#include <ranges>
#include <algorithm>


namespace tmpc
{
    using OcpVertex = unsigned int;
    using OcpEdge = unsigned int;


    class OcpTree
    {
    public:
        using vertices_size_type = unsigned int;
        using edges_size_type = unsigned int;
        using degree_size_type = unsigned int;


        explicit OcpTree(std::initializer_list<degree_size_type> branching)
        :   OcpTree(std::views::all(branching))
        {
        }
        

        /// @brief Copy constructor
        OcpTree(OcpTree const&) = default;


        /// @brief Move constructor
        OcpTree(OcpTree&&) noexcept = default;


        /// @brief Construct from branching factor sequence.
        ///
        /// @param branching a range defining branching factors for each vertex
        /// in breath-first order starting from root.
        /// The size of the range is equal to the number of vertices,
        /// and the sum of branching factors should be equal to the number of vertices minus 1.
        ///
        template <std::ranges::range Range>
        explicit OcpTree(Range const& branching)
        :   impl_ {std::make_shared<Impl>()}
        {
            impl_->numVertices_ = static_cast<vertices_size_type>(std::size(branching));

            // Check arguments
            if (impl_->numVertices_ == 0)
                TMPC_THROW_EXCEPTION(std::invalid_argument("branching must have at least 1 element"));

            // Get number of leaves
            vertices_size_type const n_leaves = std::ranges::count(
                branching.begin(), branching.end(), vertices_size_type {0});

            // Populate arrays
            impl_->firstChild_.resize(impl_->numVertices_ + 1);
            impl_->source_.resize(impl_->numVertices_ - 1);
            impl_->depth_.resize(impl_->numVertices_);
            impl_->depth_[0] = 0;
            impl_->branchVertices_.reserve(impl_->numVertices_ - n_leaves);
            impl_->leafVertices_.reserve(n_leaves);

            OcpVertex u {0}, v {1};
            
            for (auto nc : branching)
            {
                if (v + nc > impl_->numVertices_)
                    TMPC_THROW_EXCEPTION(std::invalid_argument("Inconsistent branching values"));

                impl_->firstChild_[u] = v;

                if (nc > 0)
                    impl_->branchVertices_.push_back(u);
                else
                    impl_->leafVertices_.push_back(u);

                while (nc-- > 0)
                {
                    impl_->source_[v - 1] = u;
                    impl_->depth_[v] = impl_->depth_[u] + 1;
                    ++v;
                }

                ++u;
            }

            impl_->firstChild_[u] = v;

            if (v != u)
                TMPC_THROW_EXCEPTION(std::invalid_argument("Inconsistent branching values"));

            // Init scenario count
            impl_->scenarioCount_.resize(impl_->numVertices_);

            for (auto v : vertices() | std::views::reverse)
            {
                auto const kinder = children(v);
                
                if (std::empty(kinder))
                {
                    impl_->scenarioCount_[v] = 1;
                }
                else
                {
                    impl_->scenarioCount_[v] = 0;
                    for (auto ch : kinder)
                        impl_->scenarioCount_[v] += impl_->scenarioCount_[ch];
                }
            }
        }


        /// @brief Create a linear OcpTree with num_stages vertices.
        ///    
        explicit OcpTree(size_t num_stages);


        /// @brief Create a robust MPC OcpTree.
        ///
        /// @param num_stages Number of stages of the graph. It is equal to number of control intervals + 1.
        /// \a num_stages must be positive.
        /// num_stages = 1 means one layer of nodes (consisting of the root node only);
        /// num_stages = 2 means the root node and 1 layer of nodes after it (1 time step) and so on.
        ///
        /// @param robust_horizon time step after which branching stops. robust_horizon = 0 means
        /// no branching, robust_horizon = 1 means branching only at the root node.
        ///
        /// @param branching how many children each parent node has within the robust horizon
        ///
        explicit OcpTree(size_t num_stages, size_t robust_horizon, size_t branching);


        /// @brief Forbid assignment, make OcpTree objects immutable.
        OcpTree& operator=(OcpTree const&) noexcept = delete;


        /// @brief Equality compare with other OcpTree.
        /// The OcpTree objects are equal if they describe a tree with the same structure.
        bool operator==(OcpTree const& other) const noexcept
        {
            return impl_ == other.impl_ ||
                (impl_->firstChild_.size() == other.impl_->firstChild_.size()
                && std::equal(impl_->firstChild_.begin(), impl_->firstChild_.end(), 
                    other.impl_->firstChild_.begin()));
        }


        /// @brief Number of vertices.
        vertices_size_type numVertices() const noexcept
        {
            return impl_->numVertices_;
        }


        /// @brief Number of leaves.
        vertices_size_type scenarioCount() const noexcept
        {
            return impl_->leafVertices_.size();
        }


        /// @brief Number of edges.
        edges_size_type numEdges() const noexcept
        {
            return impl_->numVertices_ - 1;
        }


        /// @brief Vertices in forward breadth-first order as a range.
        auto vertices() const noexcept
        {
            return std::views::iota(
                OcpVertex {0},
                OcpVertex {numVertices()}
            );
        }


        /// @brief Vertices with non-zero out-degree
        /// in forward breadth-first order as a range.
        auto branchVertices() const noexcept
        {
            return std::views::all(impl_->branchVertices_);
        }


        /// @brief Vertices with zero out-degree
        /// in forward breadth-first order as a range.
        auto leafVertices() const noexcept
        {
            return std::views::all(impl_->leafVertices_);
        }


        /// @brief Edges in forward breadth-first order as a range.
        auto edges() const noexcept
        {
            return std::views::iota(
                OcpEdge {0},
                OcpEdge {numEdges()}
            );
        }


        /// @brief Source vertex of a given edge.
        OcpVertex source(OcpEdge e) const
        {
            return impl_->source_.at(e);
        }


        /// @brief Target vertex of a given edge.
        OcpVertex target(OcpEdge e) const
        {
            return OcpVertex {e + 1};
        }


        /// @brief Parent edge of a given vertex.
        ///
        /// @return an std::optional<OcpEdge> containing the value equal to the
        /// parent edge of vertex \a v if \a v has parent;
        /// otherwise, an std::optional<OcpEdge> which does not contain a value.
        std::optional<OcpEdge> parentEdge(OcpVertex v) const
        {
            std::optional<OcpEdge> e;

            if (v != 0)
                e = v - 1;

            return e;
        }


        /// @brief Out-edges for a given vertex.
        ///
        /// @return a range of out-edges for vertex \a v.
        auto outEdges(OcpVertex v) const
        {
            return std::views::iota(
                OcpEdge {impl_->firstChild_.at(v) - 1}, 
                OcpEdge {impl_->firstChild_.at(v + 1) - 1}
            );
        }


        /// @brief Out-degree for a given vertex.
        ///
        /// https://en.wikipedia.org/wiki/Glossary_of_graph_theory_terms#out-degree
        ///
        /// @return out-degree (the number of out-edges) for vertex \a v.
        auto outDegree(OcpVertex v) const
        {
            return impl_->firstChild_.at(v + 1) - impl_->firstChild_.at(v);
        }


        /// @brief Number of scenarios containing given vertex.
        /// 
        /// @return number of paths from the root to a leaf node containing vertex \a v.
        auto scenarioCount(OcpVertex v) const
        {
            return impl_->scenarioCount_.at(v);
        }


        /// @brief Depth of a given vertex.
        /// 
        /// https://en.wikipedia.org/wiki/Glossary_of_graph_theory_terms#depth
        ///
        /// @return number of edges connecting vertex \a v to the root.
        auto depth(OcpVertex v) const
        {
            return impl_->depth_.at(v);
        }


        /// @brief Child vertices of a given vertex.
        /// 
        /// https://en.wikipedia.org/wiki/Glossary_of_graph_theory_terms#child
        ///
        /// @return range of vertices to which there is an out-edge from \a v.
        auto children(OcpVertex v) const
        {
            return std::views::iota(impl_->firstChild_.at(v), impl_->firstChild_.at(v + 1));
        }


    private:
        struct Impl
        {
            vertices_size_type numVertices_;

            // Maps vertex to its first child vertex.
            std::vector<OcpVertex> firstChild_;

            // Maps edge to its source vertex.
            std::vector<OcpVertex> source_;

            // Maps vertex to the number of its affected scenarios.
            std::vector<vertices_size_type> scenarioCount_;

            // Maps vertex to its depth.
            std::vector<vertices_size_type> depth_;

            // List of branch vertices.
            std::vector<OcpVertex> branchVertices_;

            // List of leaf vertices.
            std::vector<OcpVertex> leafVertices_;
        };


        std::shared_ptr<Impl> impl_;
    };

    
    inline auto num_vertices(OcpTree const& g)
    {
        return g.numVertices();
    }


    inline auto num_edges(OcpTree const& g)
    {
        return g.numEdges();
    }


    inline auto vertices(OcpTree const& g)
    {
        return g.vertices();
    }


    inline auto edges(OcpTree const& g)
    {
        return g.edges();
    }


    /// @brief Input degree of a node.
    ///
    /// Since OcpTree is always a tree, in_degree() is 0 for the root and 1 for all other nodes.
    inline OcpTree::edges_size_type in_degree(OcpVertex v, OcpTree const& g)
    {
        return g.parentEdge(v) ? 1 : 0;
    }
	
	
	inline auto out_degree(OcpVertex v, OcpTree const& g)
    {
        return g.outDegree(v);
    }


    inline OcpVertex source(OcpEdge e, OcpTree const& g)
    {
        return g.source(e);
    }


    inline OcpVertex target(OcpEdge e, OcpTree const& g)
    {
        return g.target(e);
    }
	

    /// @brief Parent of a node.
    ///
    /// Since OcpTree is always a tree, parent() is empty for the root 
    /// and equals to source(in_edges(v, g), g) for non-root nodes.
    inline std::optional<OcpVertex> parent(OcpVertex v, OcpTree const& g)
    {
        std::optional<OcpVertex> p;

        if (auto const e = g.parentEdge(v))
            p = source(*e, g);

        return p;
    }


    inline auto out_edges(OcpVertex v, OcpTree const& g)
    {
        return g.outEdges(v);
    }


    inline auto children(OcpVertex v, OcpTree const& g)
    {
        return g.children(v);
    }


    /// @brief Siblings of a vertex.
    ///
    /// Siblings of v is the set of vertices u such as parent(u, g) == parent(v, g),
    /// except v = root(g).
    /// The root is always his own only sibling.
    inline auto siblings(OcpVertex v, OcpTree const& g)
    {
        auto const p = parent(v, g);
        
        return p ? children(*p, g) : std::views::iota(v, v + 1);
    }


    /// @brief Root of the graph.
    inline OcpVertex root(OcpTree const& g)
    {
        return OcpVertex {0};
    }


    inline OcpVertex vertex(OcpTree::vertices_size_type i, OcpTree const& g)
    {
        if (i >= num_vertices(g))
            TMPC_THROW_EXCEPTION(std::out_of_range("OcpVertex index out of range"));

        return OcpVertex {i};
    }
}
