#pragma once

#include "StaticMatrix.hpp"
#include "DynamicMatrix.hpp"

#include "Eigen.hpp"

namespace tmpc
{
    template <typename MT>
    struct IdentityMatrix;

    template <typename Type, bool SO>
    struct IdentityMatrix<DynamicMatrix<Type, SO>>
    :   DynamicMatrix<Type, SO>::IdentityReturnType
    {
        typedef DynamicMatrix<Type, SO> MatrixType;
        typedef typename MatrixType::IdentityReturnType Base;

        IdentityMatrix(size_t n)
        :   Base(MatrixType::Identity(n, n))
        {            
        }
    };

    template <typename Type, size_t N, bool SO>
    struct IdentityMatrix<StaticMatrix<Type, N, N, SO>>
    :   StaticMatrix<Type, N, N, SO>::IdentityReturnType
    {
        typedef StaticMatrix<Type, N, N, SO> MatrixType;
        typedef typename MatrixType::IdentityReturnType Base;

        IdentityMatrix()
        :   Base(MatrixType::Identity())
        {            
        }
    };
}