#pragma once

#include "BinaryOp.hpp"

namespace tmpc
{
    struct Div
    {
        static std::string name()
        {
            return "Div";
        }
    
        template <typename T1, typename T2>
        static decltype(auto) eval(T1 const& a, T2 const& b)
        {
            return a / b;
        }

        template <typename T1, typename T2, typename S1, typename S2>
        static decltype(auto) diff(T1 const& a, T2 const& b, S1 const& sa, S2 const& sb)
        {
            return (sa * b - a * sb) / (b * b);
        }
    };
    

    template <typename ExprA, typename ExprB>
    decltype(auto) constexpr operator/(ExpressionBase<ExprA> const& a, ExpressionBase<ExprB> const& b)
    {
        return BinaryOp<Div, ExprA, ExprB>(a.derived(), b.derived());
    }
}