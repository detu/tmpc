#pragma once

#include <tmpc/matrix/TransposeFlag.hpp>

#include "Eigen.hpp"

#include <type_traits>

namespace tmpc :: eigen_adaptor {

template <typename VT, TransposeFlag TF>
using DenseVector = typename std::enable_if_t<
    VT::IsRowMatrix == (TF == rowVector) && IsVector<VT>::value,
    Eigen::DenseBase<MT>
>;

}