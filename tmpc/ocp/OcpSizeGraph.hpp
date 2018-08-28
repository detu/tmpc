#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/SizeT.hpp>

#include <boost/graph/adjacency_list.hpp>


namespace tmpc
{
    struct OcpVertex
    {
        OcpVertex() = default;


        OcpVertex(OcpSize const& sz)
        :   size(sz)
        {            
        }


        OcpSize size;
    };


    using OcpSizeGraph = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS, 
        OcpVertex,
        boost::property<boost::edge_index_t, size_t>
    >;


    using OcpVertexDescriptor = boost::graph_traits<OcpSizeGraph>::vertex_descriptor;
    using OcpEdgeDescriptor = boost::graph_traits<OcpSizeGraph>::edge_descriptor;


    inline auto size(OcpSizeGraph const& g)
    {
        return get(&OcpVertex::size, g);
    }


    inline std::pair<size_t, size_t> size_Q(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return std::pair(sz.nx(), sz.nx());
    }


    inline std::pair<size_t, size_t> size_R(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return std::pair(sz.nu(), sz.nu());
    }


    inline std::pair<size_t, size_t> size_S(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return std::pair(sz.nx(), sz.nu());
    }


    inline size_t size_q(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return sz.nx();
    }


    inline size_t size_r(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return sz.nu();
    }


    inline std::pair<size_t, size_t> size_C(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return std::pair(sz.nc(), sz.nx());
    }


    inline std::pair<size_t, size_t> size_D(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return std::pair(sz.nc(), sz.nu());
    }


    inline size_t size_d(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return sz.nc();
    }


    inline size_t size_x(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return sz.nx();
    }


    inline size_t size_u(OcpSizeGraph const& g, OcpVertexDescriptor v)
    {
        auto const& sz = get(size(g), v);
        return sz.nu();
    }


    inline std::pair<size_t, size_t> size_A(OcpSizeGraph const& g, OcpEdgeDescriptor e)
    {
        auto const& sz_u = get(size(g), source(e, g));
        auto const& sz_v = get(size(g), target(e, g));
        return std::pair(sz_v.nx(), sz_u.nx());
    }


    inline std::pair<size_t, size_t> size_B(OcpSizeGraph const& g, OcpEdgeDescriptor e)
    {
        auto const& sz_u = get(size(g), source(e, g));
        auto const& sz_v = get(size(g), target(e, g));
        return std::pair(sz_v.nx(), sz_u.nu());
    }


    inline size_t size_b(OcpSizeGraph const& g, OcpEdgeDescriptor e)
    {
        auto const& sz_v = get(size(g), target(e, g));
        return sz_v.nx();
    }
 

    template <typename InIterOutDegree, typename InIterOcpSize>
    inline OcpSizeGraph ocpSizeGraphFromOutDegreeList(InIterOutDegree out_degree, InIterOcpSize sz)
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
        OcpSizeGraph g(num_v);

        // Traverse the out-degree list again, create edges and set node sizes.
        od = out_degree;
        u = 0;
        v = 1;

        while (u < v)
        {
            g[u].size = *sz;

            for (size_t k = 0; k < *od; ++k)
            {
                add_edge(u, v, v - 1 /* edge index */, g);
                ++v;
            }

            ++u;
            ++od;
            ++sz;
        }
    }
}
