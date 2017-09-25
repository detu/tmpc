#pragma once

#include "UnaryOp.hpp"
#include "Identity.hpp"

#include <cmath>

namespace tmpc
{
    struct Log
    {
        static std::string name()
        {
            return "log";
        }
    
        template <typename T>
        static decltype(auto) eval(T const& arg)
        {
            return std::log(arg);
        }

        template <typename T1>
        static decltype(auto) diff(T1 const& x)
        {
            return -Identity {} / x;
        }
    };

    template <typename Expr>
    decltype(auto) constexpr log(ExpressionBase<Expr> const& x)
    {
        return UnaryOp<Log, Expr>(x.derived());
    }
}