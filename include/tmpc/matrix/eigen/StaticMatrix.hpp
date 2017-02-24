#pragma once

#include "Types.hpp"
#include "StorageOrder.hpp"
#include "MatrixAssign.hpp"
#include "StaticMatrixStorageOrder.hpp"

#include <Eigen/Dense>

namespace tmpc 
{
    template <typename Type, size_t M, size_t N, bool SO>
    using StaticMatrixBase =
        Eigen::Matrix<Type, M, N, StaticMatrixStorageOrder<M, N, SO>::value>;

    template <typename Type, size_t M, size_t N, bool SO = defaultStorageOrder>
    struct StaticMatrix
    :   StaticMatrixBase<Type, M, N, SO>
    {
        typedef StaticMatrixBase<Type, M, N, SO> Base;
        typedef Type ElementType;

        StaticMatrix()
        {            
        }

        StaticMatrix(Type const& rhs)
        {            
            this->setConstant(rhs);
        }

        StaticMatrix(initializer_list<initializer_list<Type>> list)
        {
            assign(*this, list);
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