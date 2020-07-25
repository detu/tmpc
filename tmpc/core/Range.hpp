#pragma once

#include <boost/range/irange.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/adaptor/transformed.hpp>


namespace tmpc :: views
{
    // Make it look like std::views

    template <typename Integer>
    inline auto iota(Integer first, Integer last)
    {
        return boost::irange(first, last);
    }


    template <typename Range>
    inline auto all(Range& range)
    {
        return boost::make_iterator_range(range);
    }


    // using boost::make_iterator_range;
    // using boost::adaptors::indexed;
    auto const reverse = boost::adaptors::reversed;
    // using boost::adaptors::reverse;
    auto const transform = boost::adaptors::transformed;
    // using boost::adaptors::transformed;
}