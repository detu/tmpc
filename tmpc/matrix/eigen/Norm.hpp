#pragma once

#include "Eigen.hpp"

namespace Eigen
{
    template <typename T>
    inline decltype(auto) norm(MatrixBase<T> const& m)
    {
        return m.norm();
    }

    
    template <typename T>
    inline decltype(auto) squaredNorm(MatrixBase<T> const& m)
    {
        return m.squaredNorm();
    }

    
    template <unsigned P, typename MT>
    inline decltype(auto) lpNorm(MatrixBase<MT> const& m)
    {
        return m.template lpNorm<P>();
    }

    
    template <typename MT>
    inline decltype(auto) l1Norm(MatrixBase<MT> const& m)
    {
        return m.template lpNorm<1>();
    }


    template <typename MT>
    inline decltype(auto) maxNorm(MatrixBase<MT> const& m)
    {
        return m.template lpNorm<Eigen::Infinity>();
    }
}
