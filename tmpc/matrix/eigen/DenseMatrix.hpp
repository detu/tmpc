#pragma once 

#include <tmpc/matrix/StorageOrder.hpp>

#include "Eigen.hpp"

#include <type_traits>

namespace tmpc :: eigen_adaptor {

template <typename MT, StorageOrder SO>
using DenseMatrix = typename std::enable_if_t<MT::IsRowMatrix == (SO == rowMajor), Eigen::DenseBase<MT>>;

}