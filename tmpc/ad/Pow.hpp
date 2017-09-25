#pragma once

#include "UnaryOp.hpp"

#include <cmath>

namespace tmpc
{
    struct Pow
    {
        static std::string name()
        {
            return "pow";
        }
    
        template <typename T1, typename T2>
        static decltype(auto) eval(T1 const& a, T2 const& b)
        {
            return std::pow(a, b);
        }

        template <typename T1, typename T2, typename S1, typename S2>
        static decltype(auto) diff(T1 const& a, T2 const& b, S1 const& sa, S2 const& sb)
        {
            return b * std::pow(a, b - 1.) * sa + std::log(a) * std::pow(a, b) * sb;
        }
    };

    template <typename Expr, typename T>
    decltype(auto) constexpr pow(ExpressionBase<Expr> const& a, T const& b)
    {
        return BinaryOp<Pow, Expr, T>(a.derived(), b);
    }
}