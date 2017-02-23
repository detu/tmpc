#pragma once

#include "AlignmentFlag.hpp"
#include "PaddingFlag.hpp"
#include "StorageOrder.hpp"

#include <Eigen/Dense>

namespace tmpc {

template <typename Type, bool AF, bool PF, bool SO = defaultStorageOrder>
struct CustomMatrix;

template <typename Type>
struct CustomMatrix<Type, unaligned, unpadded, rowMajor> : 
    Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
{
    typedef Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> Base;
};

template <typename Type>
struct CustomMatrix<Type, unaligned, unpadded, columnMajor> : 
    Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
{
    typedef Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> Base;
};

}