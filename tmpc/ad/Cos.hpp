#pragma once

#include "UnaryOp.hpp"

#include <cmath>

namespace tmpc
{
    struct Cos
    {
        static std::string name()
        {
            return "cos";
        }
    
        template <typename T>
        static decltype(auto) eval(T const& val)
        {
            return std::cos(val);
        }

        template <typename T1, typename T2>
        static decltype(auto) diff(T1 const& arg, T2 const& sens)
        {
            return -std::sin(arg) * sens;
        }
    };

    template <typename Expr>
    decltype(auto) constexpr cos(ExpressionBase<Expr> const& x)
    {
        return UnaryOp<Cos, Expr>(x.derived());
    }
}