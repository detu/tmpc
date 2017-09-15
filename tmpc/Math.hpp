/// \brief Some useful math functions.

#pragma once

#include <limits>

namespace tmpc
{
    template <typename T>
    inline T constexpr infinity()
    {
        return std::numeric_limits<T>::infinity();
    }
}