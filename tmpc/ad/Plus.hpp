#pragma once

#include "BinaryOp.hpp"
#include "Identity.hpp"
#include "Zero.hpp"

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

        template <typename T1, typename T2>
        static decltype(auto) diffL(T1 const& a, T2 const& b)
        {
            return Identity {};
        }

        template <typename T1, typename T2>
        static decltype(auto) diffR(T1 const& a, T2 const& b)
        {
            return Identity {};
        }
    };
    

    template <typename ExprA, typename ExprB>
    decltype(auto) constexpr operator+(ExpressionBase<ExprA> const& a, ExpressionBase<ExprB> const& b)
    {
        return BinaryOp<Plus, ExprA, ExprB>(a.derived(), b.derived());
    }


    template <typename L>
    decltype(auto) operator+(L const& a, Zero const&)
    {
        return a;
    }


    template <typename R>
    decltype(auto) operator+(Zero const&, R const& b)
    {
        return b;
    }

    
    decltype(auto) operator+(Zero const& a, Zero const& b)
    {
        return Zero {};
    }
}