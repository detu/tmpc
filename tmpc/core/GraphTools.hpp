#pragma once

#include <tmpc/core/IteratorRange.hpp>


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
}