#pragma once

#include "ForEachExpr.hpp"

namespace tmpc :: eigen_adaptor 
{
    /*------------------------------------------------------
    *
    * forEach
    *
    -------------------------------------------------------*/
    template <typename MT, typename OP>
    ForEachExpr<MT, OP> forEach(Eigen::MatrixBase<MT> const& arg, OP const& op)
    {
        return ForEachExpr<MT, OP>(arg.derived(), op);
    }
}