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

        template <typename T1, typename T2, typename S1, typename S2>
        static decltype(auto) diff(T1 const& a, T2 const& b, S1 const& sa, S2 const& sb)
        {
            return a * sb + sa * b;
        }
    };
    

    template <typename T, typename Expr>
    decltype(auto) constexpr operator*(T const& a, ExpressionBase<Expr> const& b)
    {
        return BinaryOp<Times, T, Expr>(a, b.derived());
    }
}