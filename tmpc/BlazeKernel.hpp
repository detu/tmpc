#pragma once

#include <blaze/Blaze.h>

#include <tmpc/matrix/StorageOrder.hpp>
#include <tmpc/matrix/TransposeFlag.hpp>

namespace tmpc
{
    template <typename Real_>
    struct BlazeKernel
    {
        using size_t = blaze::size_t;
        using Real = Real_;

        template <size_t M>
        using StaticVector = blaze::StaticVector<Real, M>;
        
        using DynamicVector = blaze::DynamicVector<Real>;

        template <size_t M, size_t N>
        using StaticMatrix = blaze::StaticMatrix<Real, M, N>;
        
        using DynamicMatrix = blaze::DynamicMatrix<Real>;

        template <typename T>
        using Rand = blaze::Rand<T>;
    };
}

namespace tmpc
{
    template <typename V, bool SO>
    inline decltype(auto) full(blaze::Vector<V, SO>& v)
    {
        return subvector(v, 0, size(v));
    }

    template <typename M, bool SO>
    inline decltype(auto) full(blaze::Matrix<M, SO>& m)
    {
        return submatrix(m, 0, 0, rows(m), columns(m));
    }
}