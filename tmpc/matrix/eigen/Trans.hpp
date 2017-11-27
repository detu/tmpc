#pragma once

#include "Eigen.hpp"


namespace Eigen
{
    template <typename T>
    decltype(auto) trans(MatrixBase<T> const& m)
    {
        return m.transpose();
    }
}