#pragma once

#include <tmpc/ocp/OcpSize.hpp>

#include <pair>


namespace tmpc
{
    enum class OcpQpVertexElement
    :   size_t
    {
        Q,
        R,
        S,
        q,
        r,
        lx,
        ux,
        lu,
        uu,
        C,
        D,
        ld,
        ud,
        NUM_ELEMENTS
    };


    template <OcpQpVertexElement E>
    std::pair<size_t, size_t> dimensions(OcpSize const& s);


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::Q>(OcpSize const& s)
    {
        return std::pair(s.nx(), s.nx());
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::R>(OcpSize const& s)
    {
        return std::pair(s.nu(), s.nu());
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::S>(OcpSize const& s)
    {
        return std::pair(s.nx(), s.nu());
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::q>(OcpSize const& s)
    {
        return std::pair(s.nx(), 1u);
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::r>(OcpSize const& s)
    {
        return std::pair(s.nu(), 1u);
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::lx>(OcpSize const& s)
    {
        return std::pair(s.nx(), 1u);
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::ux>(OcpSize const& s)
    {
        return std::pair(s.nx(), 1u);
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::lu>(OcpSize const& s)
    {
        return std::pair(s.nu(), 1u);
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::uu>(OcpSize const& s)
    {
        return std::pair(s.nu(), 1u);
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::C>(OcpSize const& s)
    {
        return std::pair(s.nc(), s.nx());
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::D>(OcpSize const& s)
    {
        return std::pair(s.nc(), s.nu());
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::ld>(OcpSize const& s)
    {
        return std::pair(s.nc(), 1);
    }


    template <>
    inline std::pair<size_t, size_t> dimensions<OcpQpVertexElement::ud>(OcpSize const& s)
    {
        return std::pair(s.nc(), 1);
    }


    template <OcpQpVertexElement E>
    inline size_t numElements(OcpSize const& s)
    {
        auto const dims = dimensions<E>(s);
        return dims.first * dims.second;
    }
}