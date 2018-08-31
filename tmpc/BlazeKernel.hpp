#pragma once

#include "Matrix.hpp"

#include <blaze/Blaze.h>
#include <blaze/math/IdentityMatrix.h>

namespace tmpc
{
    template <typename Real_>
    struct BlazeKernel
    {
        using size_t = blaze::size_t;
        using Real = Real_;

        template <size_t M, TransposeFlag TF>
        using StaticVector = blaze::StaticVector<Real, M, TF == rowVector ? blaze::rowVector : blaze::columnVector>;
        
        template <TransposeFlag TF>
        using DynamicVector = blaze::DynamicVector<Real, TF == rowVector ? blaze::rowVector : blaze::columnVector>;

        template <size_t M, size_t N, StorageOrder SO>
        using StaticMatrix = blaze::StaticMatrix<Real, M, N, SO == rowMajor ? blaze::rowMajor : blaze::columnMajor>;
        
        template <StorageOrder SO>
        using DynamicMatrix = blaze::DynamicMatrix<Real, SO == rowMajor ? blaze::rowMajor : blaze::columnMajor>;

        template <AlignmentFlag AF, PaddingFlag PF, StorageOrder SO>
        using CustomMatrix = blaze::CustomMatrix<
            Real,
            AF == aligned ? blaze::aligned : blaze::unaligned, 
            PF == padded ? blaze::padded : blaze::unpadded,
            SO == columnMajor ? blaze::columnMajor : blaze::rowMajor
        >;

        using IdentityMatrix = blaze::IdentityMatrix<Real>;

        template <typename MT, AlignmentFlag AF>
        using Submatrix = blaze::Submatrix<MT, AF == aligned ? blaze::aligned : blaze::unaligned>;

        template <typename VT, AlignmentFlag AF, TransposeFlag TF>
        using Subvector = blaze::Subvector<VT, AF == aligned ? blaze::aligned : blaze::unaligned, 
            TF == rowVector ? blaze::rowVector : blaze::columnVector>;

        template <typename T>
        using Rand = blaze::Rand<T>;
    };

    /*
    template <typename T>
    using KernelOf = std::enable_if_t<std::is_base_of<blaze::Matrix<T>, T>, BlazeKernel<T::ElementType>>;
    */
}

namespace blaze
{
    template <typename V, bool TF>
    inline decltype(auto) noresize(Vector<V, TF>& v)
    {
        return subvector(v, 0, size(v));
    }

    
    template <typename M, bool SO>
    inline decltype(auto) noresize(Matrix<M, SO>& m)
    {
        return submatrix(m, 0, 0, rows(m), columns(m));
    }


    template <typename V, bool TF>
    inline auto dimensions(Vector<V, TF> const& v)
    {
        return size(v);
    }


    template <typename M, bool SO>
    inline auto dimensions(Matrix<M, SO> const& m)
    {
        return std::pair(rows(m), columns(m));
    }
}