#pragma once

#include "TransposeFlag.hpp"
#include "IsVector.hpp"
#include "EigenType.hpp"

#include <type_traits>

#include <Eigen/Dense>

namespace tmpc {

template <typename VT, bool TF>
using Vector = typename std::enable_if_t<
    VT::IsRowMajor == (TF == rowVector) && IsVector<VT>::value, 
    typename EigenType<VT>::type
>;

}