#pragma once

#include "AlignmentFlag.hpp"
#include "PaddingFlag.hpp"
#include "StorageOrder.hpp"

#include <Eigen/Dense>

namespace tmpc {

template <typename Type, bool SO>
using CustomMatrixBaseUnalignedUnpadded =
    Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor>>;

template <typename Type, bool AF, bool PF, bool SO = defaultStorageOrder>
struct CustomMatrix;

template <typename Type, bool SO>
struct CustomMatrix<Type, unaligned, unpadded, SO> : 
    CustomMatrixBaseUnalignedUnpadded<Type, SO>
{
    typedef CustomMatrixBaseUnalignedUnpadded<Type, SO> Base;
    typedef Type ElementType;
};

}