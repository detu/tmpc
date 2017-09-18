#pragma once

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor
{
    template <typename T>
    decltype(auto) eval(Eigen::DenseBase<T> const& expr)
    {
        return expr.eval();
    }
}