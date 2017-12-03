#pragma once

#include <tmpc/core/ExpressionBase.hpp>


namespace tmpc
{
    ///
    /// \brief CRTP base for OCP QP expressions.
    ///
    template <typename Derived>
    class OcpQpExpressionBase
    :   public ExpressionBase<Derived>
    {
    public:
        decltype(auto) size() const
        {
            return this->derived().implSize();
        }


    protected:
		// Allow default construction and copying only as a part of a derived class;
		// otherwise, an object might be created which is not a part of Derived, 
		// and therefore calling its methods will cause undefined behavior.
		OcpQpExpressionBase() = default;
		OcpQpExpressionBase(OcpQpExpressionBase const&) = default;
    };
}