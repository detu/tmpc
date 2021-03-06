#pragma once

#include "StaticMatrix.hpp"
#include "DynamicMatrix.hpp"

#include "Eigen.hpp"

#include <tmpc/matrix/StorageOrder.hpp>

namespace tmpc :: eigen_adaptor
{
    template <typename Type, StorageOrder SO = defaultStorageOrder>
    struct IdentityMatrix
    :   DynamicMatrix<Type, SO>::IdentityReturnType
    {
        typedef DynamicMatrix<Type, SO> MatrixType;
        typedef typename MatrixType::IdentityReturnType Base;

        IdentityMatrix(size_t n)
        :   Base(MatrixType::Identity(n, n))
        {            
        }
    };
}