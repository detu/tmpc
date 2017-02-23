#pragma once

#include "StaticMatrix.hpp"
#include "DynamicMatrix.hpp"

#include <Eigen/Dense>

namespace tmpc
{
    template <typename MT>
    struct IdentityMatrix;

    template <typename Type, bool SO>
    struct IdentityMatrix<DynamicMatrix<Type, SO>>
    :   DynamicMatrix<Type, SO>::IdentityReturnType
    {
        typedef typename DynamicMatrix<Type, SO>::IdentityReturnType Base;

        IdentityMatrix(size_t n)
        :   Base(DynamicMatrix<Type, SO>::Identity(n, n))
        {            
        }
    };
}