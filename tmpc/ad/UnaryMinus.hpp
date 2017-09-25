#pragma once

#include "UnaryOp.hpp"
#include "Identity.hpp"

namespace tmpc
{
    struct UnaryMinus;

    
    decltype(auto) constexpr operator-(Identity const& I)
    {
        return UnaryOp<UnaryMinus, Identity> {I};
    }

    
    struct UnaryMinus
    {
        static std::string name()
        {
            return "unary_minus";
        }
    
        template <typename T1>
        static decltype(auto) eval(T1 const& a)
        {
            return -a;
        }

        template <typename T1>
        static decltype(auto) diff(T1 const& a)
        {
            return -Identity {};
        }
    };
    

    template <typename Expr>
    decltype(auto) constexpr operator-(ExpressionBase<Expr> const& b)
    {
        return UnaryOp<UnaryMinus, Expr> {b.derived()};
    }
}