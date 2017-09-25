#pragma once

#include "BinaryOp.hpp"
#include "Identity.hpp"
#include "Zero.hpp"


namespace tmpc
{
    struct Times
    {
        static std::string name()
        {
            return "times";
        }
    
        template <typename T1, typename T2>
        static decltype(auto) eval(T1 const& a, T2 const& b)
        {
            return a * b;
        }

        template <typename T1, typename T2>
        static decltype(auto) diffL(T1 const& a, T2 const& b)
        {
            return b;
        }

        template <typename T1, typename T2>
        static decltype(auto) diffR(T1 const& a, T2 const& b)
        {
            return a;
        }
    };
    

    template <typename L, typename R>
    decltype(auto) constexpr operator*(L const& a, ExpressionBase<R> const& b)
    {
        return BinaryOp<Times, L, R>(a, b.derived());
    }


    template <typename L>
    decltype(auto) constexpr operator*(L const& a, Identity const& b)
    {
        return a;
    }


    template <typename R>
    decltype(auto) constexpr operator*(Identity const& a, R const& b)
    {
        return b;
    }

    
    decltype(auto) constexpr operator*(Identity const& a, Identity const& b)
    {
        return Identity {};
    }


    template <typename T>
    decltype(auto) constexpr operator*(T const& a, Zero const& b)
    {
        return Zero {};
    }


    template <typename T>
    decltype(auto) constexpr operator*(Zero const& a, T const& b)
    {
        return Zero {};
    }


    decltype(auto) constexpr operator*(Zero const& a, Zero const& b)
    {
        return Zero {};
    }


    decltype(auto) constexpr operator*(Zero const&, Identity const&)
    {
        return Zero {};
    }


    decltype(auto) constexpr operator*(Identity const&, Zero const&)
    {
        return Zero {};
    }


    template <typename T>
    auto constexpr operator*(UnaryOp<UnaryMinus, Identity> const&, T const& x)
    {
        return -x;
    }
}