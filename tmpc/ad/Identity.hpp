#pragma once

#include "ExpressionBase.hpp"


namespace tmpc
{
    class Identity
    :   public ExpressionBase<Identity>
    {
    public:
        constexpr Identity()
        {            
        }

        size_t implOperationCount() const
        {
            return 0;
        }
    };
}