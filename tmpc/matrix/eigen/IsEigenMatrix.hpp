#pragma once

#include "Eigen.hpp"

#include <type_traits>

namespace tmpc :: eigen_adaptor
{
    // TODO: move to tmpc :: eigen_adaptor_ :: detail
    template <typename T>
    using IsEigenMatrix = std::is_base_of<
        Eigen::MatrixBase<
            std::remove_cv_t<T>
        >,
        std::remove_cv_t<T>
    >;
}