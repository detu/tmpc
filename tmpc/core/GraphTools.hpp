#pragma once

#include <tmpc/core/Range.hpp>


namespace tmpc
{
    template <typename Graph>
    auto verticesR(Graph const& g)
    {
        return make_iterator_range(::boost::vertices(g));
    }


    template <typename Graph>
    auto edgesR(Graph const& g)
    {
        return make_iterator_range(::boost::edges(g));
    }


    template <typename VertexDescriptor, typename Graph>
    auto adjacent_verticesR(VertexDescriptor const& v, Graph const& g)
    {
        return make_iterator_range(adjacent_vertices(v, g));
    }
}