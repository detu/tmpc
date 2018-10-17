#pragma once

#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/reversed.hpp>


namespace tmpc
{
    using boost::iterator_range;
    using boost::make_iterator_range;
    using boost::adaptors::indexed;
    using boost::adaptors::reversed;
    using boost::adaptors::reverse;
}