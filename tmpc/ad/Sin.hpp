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
        static decltype(auto) eval(T const& arg)
        {
            return std::sin(arg);
        }

        template <typename T1>
        static decltype(auto) diff(T1 const& x)
        {
            return cos(x);
        }
    };

    template <typename Expr>
    decltype(auto) constexpr sin(ExpressionBase<Expr> const& x)
    {
        return UnaryOp<Sin, Expr>(x.derived());
    }
}