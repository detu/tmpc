#pragma once

#include "UnaryOp.hpp"

#include <cmath>

namespace tmpc
{
    struct Sin
    {
        static std::string name()
        {
            return "sin";
        }
    
        template <typename T>
        static decltype(auto) eval(T const& val)
        {
            return std::sin(val);
        }
    };

    template <typename Expr>
    decltype(auto) constexpr sin(ExpressionBase<Expr> const& x)
    {
        return UnaryOp<Sin, Expr>(x.derived());
    }
}