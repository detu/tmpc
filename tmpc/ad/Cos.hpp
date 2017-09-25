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
    };

    template <typename Expr>
    decltype(auto) constexpr cos(ExpressionBase<Expr> const& x)
    {
        return UnaryOp<Cos, Expr>(x.derived());
    }
}