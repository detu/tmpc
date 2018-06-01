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
            auto const size = get(OcpSizeProperty_t(), graph);

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

                size_vec.push_back(get(size, *v));
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
            size_ = get(OcpSizeProperty_t(), graph_);
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


    private:
        using SizeMap = typename boost::property_map<OcpSizeGraph, OcpSizeProperty_t>::type;
        using VertexIdMap = boost::property_map<OcpSizeGraph, boost::vertex_index_t>::type;
        using EdgeIdMap = boost::property_map<OcpSizeGraph, boost::edge_index_t>::type;
        using ProblemVertexMap = boost::iterator_property_map<typename QpWorkspace::ProblemIterator, VertexIdMap>;
        using ProblemEdgeMap = boost::iterator_property_map<typename QpWorkspace::ProblemIterator, EdgeIdMap>;

        OcpSizeGraph graph_;
        SizeMap size_;
        VertexIdMap vertexId_;
        EdgeIdMap edgeId_;

        QpWorkspace workspace_;
        ProblemVertexMap problemVertex_;
        ProblemEdgeMap problemEdge_;
    };
}