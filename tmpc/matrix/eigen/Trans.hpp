#pragma once

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor
{
    template <typename T>
    decltype(auto) trans(Eigen::MatrixBase<T> const& m)
    {
        return m.transpose();
    }
}