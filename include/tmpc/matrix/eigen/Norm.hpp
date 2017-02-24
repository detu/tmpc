#pragma once

#include <Eigen/Dense>

namespace tmpc
{
    template <typename T>
    decltype(auto) norm(Eigen::MatrixBase<T> const& m)
    {
        return m.norm();
    }

    template <typename T>
    decltype(auto) squaredNorm(Eigen::MatrixBase<T> const& m)
    {
        return m.squaredNorm();
    }
}