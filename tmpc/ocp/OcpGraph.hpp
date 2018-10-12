#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/core/Graph.hpp>

#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/compressed_sparse_row_graph.hpp>


namespace tmpc
{
    //using OcpGraph = boost::compressed_sparse_row_graph<boost::bidirectionalS>;

    using OcpGraph = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::bidirectionalS,
        boost::no_property, boost::property<boost::edge_index_t, size_t>
    >;


    using OcpVertexDescriptor = boost::graph_traits<OcpGraph>::vertex_descriptor;
    using OcpEdgeDescriptor = boost::graph_traits<OcpGraph>::edge_descriptor;


    /// Create a tree-structured OcpGraph from a list defined by the out_gedree iterator,
    /// which is the number out-edges of each node, in breadth-first order.
    ///
    template <typename InIterOutDegree>
    inline OcpGraph ocpGraphFromOutDegree(InIterOutDegree out_degree)
    {
        // Traverse the out-degree list and calculate the total number of vertices.
        auto od = out_degree;
        size_t u = 0;
        size_t v = 1;

        while (u < v)
        {
            v += *od;
            ++u;
            ++od;
        }

        size_t const num_v = u;

        // Create the graph.
        OcpGraph g(num_v);

        // Traverse the out-degree list again, create edges and set node sizes.
        od = out_degree;
        u = 0;
        v = 1;

        while (u < v)
        {
            for (size_t k = 0; k < *od; ++k)
            {
                add_edge(u, v, v - 1 /* edge index */, g);
                ++v;
            }

            ++u;
            ++od;
        }

        return g;
    }


    /// Create a linear OcpGraph with n nodes.
    ///
    OcpGraph ocpGraphLinear(size_t n);


    /// Vertex index property map.
    inline auto vertexIndex(OcpGraph const& g)
    {
        return get(boost::vertex_index, g);
    }
}
