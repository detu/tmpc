#pragma once

#include "BinaryOp.hpp"

namespace tmpc
{
    struct Minus
    {
        static std::string name()
        {
            return "minus";
        }
    
        template <typename T1, typename T2>
        static decltype(auto) eval(T1 const& a, T2 const& b)
        {
            return a - b;
        }

        template <typename T1, typename T2>
        static decltype(auto) diffL(T1 const& a, T2 const& b)
        {
            return Identity {};
        }

        template <typename T1, typename T2>
        static decltype(auto) diffR(T1 const& a, T2 const& b)
        {
            return -Identity {};
        }
    };
    

    template <typename T, typename Expr>
    decltype(auto) constexpr operator-(T const& a, ExpressionBase<Expr> const& b)
    {
        return BinaryOp<Minus, T, Expr>(a, b.derived());
    }
}