#pragma once

#include "StaticMatrix.hpp"
#include "DynamicMatrix.hpp"

#include "Eigen.hpp"

namespace tmpc
{
    template <typename Type, bool SO = defaultStorageOrder>
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