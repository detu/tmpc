#include "StorageOrder.hpp"
#include "EigenType.hpp"

#include <Eigen/Dense>

#include <type_traits>

namespace tmpc {

template <typename MT, bool SO>
using Matrix = std::enable_if_t<
    MT::IsRowMajor == (SO == rowMajor), 
    Eigen::MatrixBase<typename EigenType<MT>::type>
>;

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