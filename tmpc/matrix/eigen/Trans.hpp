#pragma once

#include "Eigen.hpp"

namespace tmpc
{
    template <typename T>
    decltype(auto) trans(Eigen::MatrixBase<T> const& m)
    {
        return m.transpose();
    }
}