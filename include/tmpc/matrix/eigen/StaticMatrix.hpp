#pragma once

#include "Types.hpp"
#include "StorageOrder.hpp"
#include "MatrixAssign.hpp"
#include "StaticMatrixStorageOrder.hpp"

#include <Eigen/Dense>

namespace tmpc 
{

    template <typename Type, size_t M, size_t N, bool SO = defaultStorageOrder>
    struct StaticMatrix
    :   Eigen::Matrix<Type, M, N, StaticMatrixStorageOrder<M, N, SO>::value>
    {
        typedef Eigen::Matrix<Type, M, N, StaticMatrixStorageOrder<M, N, SO>::value> Base;

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

}