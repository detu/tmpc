#pragma once

#include <Eigen/Dense>

namespace tmpc
{
    template <typename T>
    decltype(auto) trans(Eigen::MatrixBase<T> const& m)
    {
        return m.transpose();
    }
}