#include "StorageOrder.hpp"

#include <Eigen/Dense>

#include <type_traits>

namespace tmpc {

template <typename MT, bool SO>
using Matrix = typename std::enable_if<MT::IsRowMatrix == (SO == rowMajor), Eigen::MatrixBase<MT>>::type;

template <typename MT>
size_t rows(Eigen::MatrixBase<MT> const& m)
{
    return m.rows();
}

template <typename MT>
size_t columns(Eigen::MatrixBase<MT> const& m)
{
    return m.cols();
}

}