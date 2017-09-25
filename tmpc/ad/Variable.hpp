#pragma once

#include "ExpressionBase.hpp"

#include <tuple>

namespace tmpc
{
    template <std::size_t N_>
    class Variable
    :   public ExpressionBase<Variable<N_>>
    {
    public:
        static auto constexpr N = N_;
    };
    
    
    template <std::size_t N, class... Types>
    decltype(auto) eval(Variable<N> const&, std::tuple<Types...> const& arg)
    {
        return std::get<N>(arg);
    }
}