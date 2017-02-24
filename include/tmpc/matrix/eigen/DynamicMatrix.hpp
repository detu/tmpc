#pragma once 

#include "StorageOrder.hpp"
#include "EigenType.hpp"
#include "InitializerList.hpp"
#include "MatrixAssign.hpp"

#include <Eigen/Dense>

namespace tmpc {

template <typename Type, bool SO>
using DynamicMatrixBase =
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor>;

template <typename Type, bool SO = defaultStorageOrder>
struct DynamicMatrix
:   DynamicMatrixBase<Type, SO>
{
    typedef DynamicMatrixBase<Type, SO> Base;
    typedef Type ElementType;

    DynamicMatrix(size_t M, size_t N)
    :   Base(M, N)
    {        
    }

    DynamicMatrix(size_t M, size_t N, Type const& val)
    :   Base(M, N)
    {
        this->setConstant(val);
    }

    DynamicMatrix(initializer_list<initializer_list<Type>> list)
    :   Base(list.size(), determineColumns(list))
    {
        assign(*this, list);
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