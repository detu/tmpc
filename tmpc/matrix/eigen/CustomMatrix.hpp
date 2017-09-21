#pragma once

#include "PaddingFlag.hpp"
#include <tmpc/matrix/StorageOrder.hpp>
#include <tmpc/matrix/AlignmentFlag.hpp>

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor {

template <typename Type, StorageOrder SO>
using CustomMatrixBaseUnalignedUnpadded =
    Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor>>;

template <typename Type, AlignmentFlag AF, bool PF, StorageOrder SO = defaultStorageOrder>
struct CustomMatrix;

template <typename Type, StorageOrder SO>
struct CustomMatrix<Type, unaligned, unpadded, SO> : 
    CustomMatrixBaseUnalignedUnpadded<Type, SO>
{
    typedef CustomMatrixBaseUnalignedUnpadded<Type, SO> Base;
    typedef Type ElementType;
};

}