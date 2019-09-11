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


        template <AlignmentFlag AF, PaddingFlag PF, TransposeFlag TF>
        using CustomVector = blaze::CustomVector<
            Real,
            AF == aligned ? blaze::aligned : blaze::unaligned, 
            PF == padded ? blaze::padded : blaze::unpadded,
            TF == columnVector ? blaze::columnVector : blaze::rowVector
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
}
