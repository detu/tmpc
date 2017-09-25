#pragma once

#include <tuple>

namespace tmpc
{
    template <typename T, class... Types>
    decltype(auto) constexpr eval(T const& x, std::tuple<Types...> const& arg)
    {
        return x;
    }
}