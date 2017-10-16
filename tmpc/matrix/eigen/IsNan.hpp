#pragma once

#include "Eigen.hpp"

namespace Eigen
{
    // Returns true if at least one element of a Matrix/Vector is a NaN.
    template <typename T>
    bool isnan(MatrixBase<T> const& expr)
    {
        return isnan(expr.array()).any();
    }
}