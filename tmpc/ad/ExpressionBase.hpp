#pragma once

#include <tmpc/SizeT.hpp>

namespace tmpc
{
    template <typename Derived>
    class ExpressionBase
    {
    public:
        constexpr Derived& derived()
        {
            return static_cast<Derived&>(*this);
        }
    
        constexpr Derived const& derived() const
        {
            return static_cast<Derived const&>(*this);
        }

        constexpr size_t operationCount() const
        {
            return derived().implOperationCount();
        }

    protected:
        constexpr ExpressionBase() = default;

        ExpressionBase(ExpressionBase const&)
        {            
        }

        ExpressionBase(ExpressionBase &&)
        {            
        }
    };


    template <typename Expr>
    constexpr size_t operationCount(ExpressionBase<Expr> const& expr)
    {
        return expr.operationCount();
    }


    constexpr size_t operationCount(double)
    {
        return 0;
    }
}