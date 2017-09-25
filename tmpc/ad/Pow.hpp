#pragma once

#include "UnaryOp.hpp"
#include "Log.hpp"

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

        template <typename T1, typename T2>
        static decltype(auto) diffL(T1 const& a, T2 const& b)
        {
            return b * pow(a, b - 1);
        }

        template <typename T1, typename T2>
        static decltype(auto) diffR(T1 const& a, T2 const& b)
        {
            return log(a) * pow(a, b);
        }
    };

    template <typename Expr, typename T>
    decltype(auto) constexpr pow(ExpressionBase<Expr> const& a, T const& b)
    {
        return BinaryOp<Pow, Expr, T>(a.derived(), b);
    }
}