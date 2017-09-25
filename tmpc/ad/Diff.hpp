#pragma once

#include "Forward.hpp"
#include "Zero.hpp"


namespace tmpc
{
    template <typename T, std::size_t N>
    decltype(auto) constexpr diff(T const& x, Variable<N> const& var)
    {
        return Zero {};
    }
}