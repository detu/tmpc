#pragma once

#include <tmpc/ocp/OcpSizeGraph.hpp>
#include <tmpc/core/IteratorRange.hpp>

#include <boost/graph/copy.hpp>


namespace tmpc
{
    namespace detail
    {
        template <typename Graph>
        inline std::vector<OcpSize> sizeGraphToSizeVector(Graph const& graph)
        {
            auto const n_stages = num_vertices(graph);
            auto const vert = make_iterator_range(vertices(graph));

            std::vector<OcpSize> size_vec;
            size_vec.reserve(n_stages);

            size_t vertex_index = 0;
            for (auto v = vert.begin(); v != vert.end(); ++v, ++vertex_index)
            {
                auto adj_vert = make_iterator_range(adjacent_vertices(*v, graph));

                if (adj_vert.size() > 1)
                    throw std::invalid_argument("TreeQpWorkspaceAdaptor does not support tree QP structures");

                if (adj_vert.size() == 0 && std::next(v) != vert.end())
                    throw std::invalid_argument("TreeQpWorkspaceAdaptor does not support disconnected graphs");
                
                if (adj_vert.size() == 1 && adj_vert.front() != vertex_index + 1)
                    throw std::invalid_argument("TreeQpWorkspaceAdaptor: graph nodes must be sequentially connected");

                size_vec.push_back(graph[*v].size);
            }

            return size_vec;
        }
    }


    template <typename QpWorkspace>
    class TreeQpWorkspaceAdaptor
    {
    public:
        using Kernel = typename QpWorkspace::Kernel;
        using Real = typename QpWorkspace::Real;

        /**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <typename Graph>
		explicit TreeQpWorkspaceAdaptor(Graph const& graph)
        :   workspace_ {detail::sizeGraphToSizeVector(graph)}
		{
			copy_graph(graph, graph_);
            size_ = get(&OcpVertex::size, std::as_const(graph_));
            vertexId_ = get(boost::vertex_index, graph_);
            edgeId_ = get(boost::edge_index, graph_);
            
            problemVertex_ = ProblemVertexMap(workspace_.problem().begin(), vertexId_);            
            problemEdge_ = ProblemEdgeMap(workspace_.problem().begin(), edgeId_);
		}


        auto const& graph() const
        {
            return graph_;
        }


        auto const& size() const
        {
            return size_;
        }


        auto const& problemVertex() const
        {
            return problemVertex_;
        }


        auto const& problemEdge() const
        {
            return problemEdge_;
        }


        auto const& vertexIndex() const
        {
            return vertexId_;
        }


        auto const& edgeIndex() const
        {
            return edgeId_;
        }


        void solve()
        {
            workspace_.solve();
        }


    private:
        class ProblemVertexIterator
		:	public boost::iterator_adaptor<
				ProblemVertexIterator	// derived
			,	typename QpWorkspace::ProblemIterator	// base
			,	OcpQpVertexBase<typename QpWorkspace::Stage>&	// value
			,	boost::random_access_traversal_tag	// category of traversal
			,	OcpQpVertexBase<typename QpWorkspace::Stage>&	// reference
			>
		{
		public:
			ProblemVertexIterator() = default;
			
			ProblemVertexIterator(typename ProblemVertexIterator::base_type const& p)
			:	ProblemVertexIterator::iterator_adaptor_(p)
			{			
			}
		};


        class ProblemEdgeIterator
		:	public boost::iterator_adaptor<
				ProblemEdgeIterator	// derived
			,	typename QpWorkspace::ProblemIterator	// base
			,	OcpQpEdgeBase<typename QpWorkspace::Stage>&	// value
			,	boost::random_access_traversal_tag	// category of traversal
			,	OcpQpEdgeBase<typename QpWorkspace::Stage>&	// reference
			>
		{
		public:
			ProblemEdgeIterator() = default;
			
			ProblemEdgeIterator(typename ProblemEdgeIterator::base_type const& p)
			:	ProblemEdgeIterator::iterator_adaptor_(p)
			{			
			}
		};


        using SizeMap = typename boost::property_map<OcpSizeGraph, OcpSize OcpVertex::*>::const_type;
        using VertexIdMap = boost::property_map<OcpSizeGraph, boost::vertex_index_t>::type;
        using EdgeIdMap = boost::property_map<OcpSizeGraph, boost::edge_index_t>::type;
        using ProblemVertexMap = boost::iterator_property_map<ProblemVertexIterator, VertexIdMap>;
        using ProblemEdgeMap = boost::iterator_property_map<ProblemEdgeIterator, EdgeIdMap>;


        OcpSizeGraph graph_;
        SizeMap size_;
        VertexIdMap vertexId_;
        EdgeIdMap edgeId_;

        QpWorkspace workspace_;
        ProblemVertexMap problemVertex_;
        ProblemEdgeMap problemEdge_;
    };
}