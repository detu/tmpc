#pragma once

#include <tmpc/matrix/EigenAdaptor.hpp>

namespace tmpc
{
    template <typename Real_>
    struct EigenKernel
    {
        using size_t = std::size_t;
        using Real = Real_;

        template <size_t M>
        using StaticVector = eigen_adaptor::StaticVector<Real, M>;
        
        using DynamicVector = eigen_adaptor::DynamicVector<Real>;

        template <size_t M, size_t N>
        using StaticMatrix = eigen_adaptor::StaticMatrix<Real, M, N>;
        
        using DynamicMatrix = eigen_adaptor::DynamicMatrix<Real>;
    };
}
