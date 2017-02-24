#pragma once

#include "Types.hpp"
#include "StorageOrder.hpp"
#include "MatrixAssign.hpp"
#include "Matrix.hpp"

#include <Eigen/Dense>

namespace tmpc 
{
    template <typename Type, size_t M, size_t N, bool SO = defaultStorageOrder>
    struct StaticMatrix
    :   Matrix<StaticMatrix<Type, M, N, SO>, SO>
    {
        typedef Matrix<StaticMatrix<Type, M, N, SO>, SO> Base;

        StaticMatrix()
        {            
        }

        StaticMatrix(Type const& rhs)
        :   Base(rhs)
        {            
        }

        StaticMatrix(initializer_list<initializer_list<Type>> list)
        :   Base(list)
        {
        }

        template <typename T>
        StaticMatrix(Eigen::MatrixBase<T> const& rhs)
        :   Base(rhs)
        {        
        }
    };

    template <typename Type, size_t M, size_t N, bool SO>
    struct EigenType<StaticMatrix<Type, M, N, SO>>
    {
        typedef typename StaticMatrix<Type, M, N, SO>::Base type;
    };
}