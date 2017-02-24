#pragma once 

#include "StorageOrder.hpp"

#include "Eigen.hpp"

#include <type_traits>

namespace tmpc {

template <typename MT, bool SO>
using DenseMatrix = typename std::enable_if_t<MT::IsRowMatrix == (SO == rowMajor), Eigen::DenseBase<MT>>;

}