#pragma once

#include "ExpressionBase.hpp"
#include "Eval.hpp"
#include "Diff.hpp"
#include "Zero.hpp"
#include "Identity.hpp"


namespace tmpc
{
    template <std::size_t N_>
    class Variable final
    :   public ExpressionBase<Variable<N_>>
    {
    public:
        constexpr Variable()
        {
        }

        constexpr Variable(Variable const&)
        {
        }

        size_t implOperationCount() const
        {
            return 0;
        }

        static auto constexpr N = N_;
    };
    
    
    template <std::size_t N, class... Types>
    decltype(auto) eval(Variable<N> const&, std::tuple<Types...> const& arg)
    {
        return std::get<N>(arg);
    }


    template <std::size_t N1, std::size_t N2>
    decltype(auto) constexpr diff(Variable<N1> const& v, Variable<N2> const& u)
    {
        return Zero {};
    }


    template <std::size_t N>
    decltype(auto) constexpr diff(Variable<N> const& v, Variable<N> const& u)
    {
        return Identity {};
    }
}