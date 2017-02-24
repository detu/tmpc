#pragma once

#include <Eigen/Dense>

namespace tmpc
{
    template <typename T>
    decltype(auto) eval(Eigen::DenseBase<T> const& expr)
    {
        return expr.eval();
    }
}