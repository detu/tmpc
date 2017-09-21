#pragma once

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor
{
   
}

namespace Eigen
{
    template <typename T>
    decltype(auto) trans(MatrixBase<T> const& m)
    {
        return m.transpose();
    }
}