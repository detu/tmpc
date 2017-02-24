#pragma once 

#include "StorageOrder.hpp"
#include "EigenType.hpp"
#include "InitializerList.hpp"
#include "MatrixAssign.hpp"
#include "Matrix.hpp"

#include <Eigen/Dense>

namespace tmpc 
{
    template <typename Type, bool SO = defaultStorageOrder>
    struct DynamicMatrix
    :   Matrix<DynamicMatrix<Type, SO>, SO>
    {
        typedef Matrix<DynamicMatrix<Type, SO>, SO> Base;

        DynamicMatrix(size_t M, size_t N)
        :   Base(M, N)
        {        
        }

        DynamicMatrix(size_t M, size_t N, Type const& val)
        :   Base(M, N, val)
        {
        }

        DynamicMatrix(initializer_list<initializer_list<Type>> list)
        :   Base(list)
        {
        }

        template <typename T>
        DynamicMatrix(Eigen::MatrixBase<T> const& rhs)
        :   Base(rhs)
        {        
        }
    };

    template <typename Type, bool SO>
    struct EigenType<DynamicMatrix<Type, SO>>
    {
        typedef typename DynamicMatrix<Type, SO>::Base type;
    };

}