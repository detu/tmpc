#pragma once

#include <Eigen/Dense>

#include <type_traits>

namespace tmpc
{
    template <typename T>
    using IsEigenMatrix = std::is_base_of<Eigen::MatrixBase<T>, T>;
}