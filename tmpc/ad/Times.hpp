#pragma once

#include "BinaryOp.hpp"

namespace tmpc
{
    struct Times
    {
        static std::string name()
        {
            return "times";
        }
    
        template <typename T1, typename T2>
        static decltype(auto) eval(T1 const& a, T2 const& b)
        {
            return a * b;
        }
    };
    

    template <typename T, typename Expr>
    decltype(auto) constexpr operator*(T const& a, ExpressionBase<Expr> const& b)
    {
        return BinaryOp<Times, T, Expr>(a, b.derived());
    }
}