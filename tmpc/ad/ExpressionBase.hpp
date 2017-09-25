#pragma once

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
    };
}