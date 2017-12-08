#pragma once

#include "IteratorRange.hpp"
#include "TransformIterator.hpp"


namespace tmpc
{
    template <typename Iterator, typename UnaryFunction>
    inline decltype(auto) make_transform_iterator_range(Iterator first, Iterator last, UnaryFunction f)
    {
        return make_iterator_range(make_transform_iterator(first, f), make_transform_iterator(last, f));
    }
}