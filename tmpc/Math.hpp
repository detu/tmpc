/// \brief Some useful math functions.

#pragma once

#include <limits>

namespace tmpc
{
    template <typename T>
    inline T constexpr inf()
    {
        return std::numeric_limits<T>::infinity();
    }

    template <typename T>
    inline T constexpr sNaN()
    {
        return std::numeric_limits<T>::signaling_NaN();
    }
}