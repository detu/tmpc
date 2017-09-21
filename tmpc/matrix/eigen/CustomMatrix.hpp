#pragma once

#include <tmpc/matrix/PaddingFlag.hpp>
#include <tmpc/matrix/StorageOrder.hpp>
#include <tmpc/matrix/AlignmentFlag.hpp>

#include "Eigen.hpp"
#include "Types.hpp"

namespace tmpc :: eigen_adaptor {

template <typename Type, StorageOrder SO>
using CustomMatrixBaseUnalignedUnpadded =
    Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor>>;

template <typename Type, AlignmentFlag AF, PaddingFlag PF, StorageOrder SO = defaultStorageOrder>
struct CustomMatrix;

template <typename Type, StorageOrder SO>
struct CustomMatrix<Type, unaligned, unpadded, SO> : 
    CustomMatrixBaseUnalignedUnpadded<Type, SO>
{
    typedef CustomMatrixBaseUnalignedUnpadded<Type, SO> Base;
    typedef Type ElementType;

    CustomMatrix(ElementType * ptr, size_t m, size_t n)
    :   Base {ptr, static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(n)}
    {        
    }

    CustomMatrix& operator=(ElementType const& rhs)
    {
        this->setConstant(rhs);
        return *this;
    }
};

}