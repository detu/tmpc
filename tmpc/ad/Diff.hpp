#pragma once

#include "Forward.hpp"


namespace tmpc
{
    template <typename T, std::size_t N, typename S, class... Types>
    decltype(auto) constexpr diff(T const& x, Variable<N> const& var, std::tuple<Types...> const& arg, S const& sens)
    {
        return S {};
    }
}