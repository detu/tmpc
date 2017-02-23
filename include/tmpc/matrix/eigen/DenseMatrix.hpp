#pragma once 

#include "StorageOrder.hpp"

#include <Eigen/Dense>

#include <type_traits>

namespace tmpc {

template <typename MT, bool SO>
using DenseMatrix = typename std::enable_if<MT::IsRowMatrix == (SO == rowMajor), Eigen::DenseBase<MT>>::type;

}