#pragma once

#include "TransposeFlag.hpp"

#include <Eigen/Dense>

#include <type_traits>

namespace tmpc {

template <typename VT, bool TF>
using DenseVector = typename std::enable_if_t<
    VT::IsRowMatrix == (TF == rowVector) && IsVector<VT>::value,
    Eigen::DenseBase<MT>
>;

}