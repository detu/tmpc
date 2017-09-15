#pragma once

#include "Eigen.hpp"

namespace Eigen
{
    template <typename T>
    decltype(auto) norm(MatrixBase<T> const& m)
    {
        return m.norm();
    }

    template <typename T>
    decltype(auto) squaredNorm(MatrixBase<T> const& m)
    {
        return m.squaredNorm();
    }

    template <unsigned P, typename MT>
    decltype(auto) lpNorm(MatrixBase<MT> const& m)
    {
        return m.template lpNorm<P>();
    }

    template <typename MT>
    decltype(auto) l1Norm(MatrixBase<MT> const& m)
    {
        return m.template lpNorm<1>();
    }
}