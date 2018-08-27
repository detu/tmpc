#pragma once

#include <tmpc/ocp/OcpSize.hpp>

#include <pair>


namespace tmpc
{
    enum class OcpQpEdgeElement
    :   size_t
    {
        A,
        B,
        b
    };


    template <OcpQpEdgeElement E>
    std::pair<size_t, size_t> dimensions(OcpSize const& s);


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpEdgeElement::A>(OcpSizeGraph const& g, OcpEdgeDescriptor e)
    {
        auto const& sz_u = get(size(g), source(e, g));
        auto const& sz_v = get(size(g), target(e, g));
        return std::pair(sz_v.nx(), sz_u.nx());
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpEdgeElement::B>(OcpSizeGraph const& g, OcpEdgeDescriptor e)
    {
        auto const& sz_u = get(size(g), source(e, g));
        auto const& sz_v = get(size(g), target(e, g));
        return std::pair(sz_v.nx(), sz_u.nu());
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpEdgeElement::b>(OcpSizeGraph const& g, OcpEdgeDescriptor e)
    {
        auto const& sz_v = get(size(g), target(e, g));
        return std::pair(sz_v.nx(), 1);
    }


    template <OcpQpEdgeElement E>
    inline size_t numElements(OcpSizeGraph const& g, OcpEdgeDescriptor e)
    {
        auto const dims = dimensions<E>(g, e);
        return dims.first * dims.second;
    }
}