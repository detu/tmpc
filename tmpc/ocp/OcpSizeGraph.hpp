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
}
