#pragma once

#include "BinaryOp.hpp"

namespace tmpc
{
    struct Plus
    {
        static std::string name()
        {
            return "plus";
        }
    
        template <typename T1, typename T2>
        static decltype(auto) eval(T1 const& a, T2 const& b)
        {
            return a + b;
        }
    };
    

    template <typename ExprA, typename ExprB>
    decltype(auto) constexpr operator+(ExpressionBase<ExprA> const& a, ExpressionBase<ExprB> const& b)
    {
        return BinaryOp<Plus, ExprA, ExprB>(a.derived(), b.derived());
    }
}