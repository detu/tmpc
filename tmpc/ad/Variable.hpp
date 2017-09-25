#pragma once

#include "ExpressionBase.hpp"
#include "Eval.hpp"
#include "Diff.hpp"


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


    template <std::size_t N1, std::size_t N2, typename S, class... Types>
    decltype(auto) constexpr diff(Variable<N1> const& v, Variable<N2> const& u, std::tuple<Types...> const& arg, S const& sens)
    {
        return S {};
    }


    template <std::size_t N, typename S, class... Types>
    decltype(auto) constexpr diff(Variable<N> const& v, Variable<N> const& u, std::tuple<Types...> const& arg, S const& sens)
    {
        return sens;
    }
}