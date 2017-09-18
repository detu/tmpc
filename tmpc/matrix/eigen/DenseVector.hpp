#pragma once

#include "TransposeFlag.hpp"

#include "Eigen.hpp"

#include <type_traits>

namespace tmpc :: eigen_adaptor {

template <typename VT, bool TF>
using DenseVector = typename std::enable_if_t<
    VT::IsRowMatrix == (TF == rowVector) && IsVector<VT>::value,
    Eigen::DenseBase<MT>
>;

}