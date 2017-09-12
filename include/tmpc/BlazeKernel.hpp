#pragma once

#include <blaze/Blaze.h>

namespace tmpc
{
    template <typename Scalar_>
    struct BlazeKernel
    {
        using size_t = blaze::size_t;
        using Scalar = Scalar_;

        template <size_t M>
        using StaticVector = blaze::StaticVector<Scalar, M>;
        
        using DynamicVector = blaze::DynamicVector<Scalar>;

        template <size_t M, size_t N>
        using StaticMatrix = blaze::StaticMatrix<Scalar, M, N>;
        
        using DynamicMatrix = blaze::DynamicMatrix<Scalar>;
    };
}

namespace blaze
{
    template <typename V, bool SO>
    inline decltype(auto) full(Vector<V, SO>& v)
    {
        return subvector(v, 0, size(v));
    }

    template <typename M, bool SO>
    inline decltype(auto) full(Matrix<M, SO>& m)
    {
        return submatrix(m, 0, 0, rows(m), columns(m));
    }
}