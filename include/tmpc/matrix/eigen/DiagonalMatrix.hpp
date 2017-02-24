#pragma once

#include <Eigen/Dense>

#include <type_traits>

namespace tmpc {

/*------------------------------------------------------
 *
 * DiagonalMatrix adaptor
 *
 -------------------------------------------------------*/
template <typename MT>
using DiagonalMatrix = typename std::enable_if<MT::RowsAtCompileTime == MT::ColsAtCompileTime, 
    Eigen::DiagonalMatrix<typename MT::Scalar, MT::RowsAtCompileTime>>::type;

}