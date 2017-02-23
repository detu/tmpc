#pragma once

#include "TransposeFlag.hpp"
#include "IsVector.hpp"

#include <type_traits>

#include <Eigen/Dense>

namespace tmpc {

template <typename VT, bool TF>
using Vector = typename std::enable_if<
    VT::IsRowMatrix == (TF == rowVector) && IsVector<VT>::value, 
    Eigen::MatrixBase<VT>
>::type;

}