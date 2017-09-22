#pragma once

#include "Eigen.hpp"

namespace Eigen
{
    template <typename T>
    decltype(auto) eval(DenseBase<T> const& expr)
    {
        return expr.eval();
    }
}