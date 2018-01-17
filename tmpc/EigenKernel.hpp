#pragma once

#include <tmpc/matrix/EigenAdaptor.hpp>
#include <tmpc/Matrix.hpp>

#include <type_traits>

namespace tmpc
{
    template <typename Real_>
    struct EigenKernel
    {
        using size_t = eigen_adaptor::size_t;
        using Real = Real_;

        template <size_t M, TransposeFlag TF>
        using StaticVector = eigen_adaptor::StaticVector<Real, M, TF>;
        
        template <TransposeFlag TF>
        using DynamicVector = eigen_adaptor::DynamicVector<Real, TF>;

        template <size_t M, size_t N, StorageOrder SO>
        using StaticMatrix = eigen_adaptor::StaticMatrix<Real, M, N, SO>;
        
        template <StorageOrder SO>
        using DynamicMatrix = eigen_adaptor::DynamicMatrix<Real, SO>;

        template <AlignmentFlag AF, PaddingFlag PF, StorageOrder SO>
        using CustomMatrix = eigen_adaptor::CustomMatrix<Real, AF, PF, SO>;

        using IdentityMatrix = eigen_adaptor::IdentityMatrix<Real>;

        template <typename MT, AlignmentFlag AF>
        using Submatrix = eigen_adaptor::Submatrix<MT, AF>;

        template <typename VT, AlignmentFlag AF, TransposeFlag TF>
        using Subvector = eigen_adaptor::Subvector<VT, AF, TF>;

        template <typename T>
        using Rand = eigen_adaptor::Rand<T>;
    };

    template <typename T>
    using KernelOf = std::enable_if_t<std::is_base_of<Eigen::EigenBase<T>, T>::value, EigenKernel<typename T::Scalar>>;
}


namespace Eigen
{
    template <typename MT>
    inline void randomize(MatrixBase<MT>& m)
    {
        m = MT::Random(m.rows(), m.cols());
    }
}