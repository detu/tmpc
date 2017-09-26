#pragma once

#include "Eigen.hpp"

#include <type_traits>

namespace tmpc :: eigen_adaptor 
{
    /*
    template <typename T, typename Enable = void>
    struct IsVector
    {
        static bool constexpr value = false;
    };

    template <typename T>
    struct IsVector<T, std::enable_if_t<std::is_base_of<Eigen::DenseBase<T>, typename T::Base>::value>>
    {
        static bool constexpr value = T::RowsAtCompileTime == 1 || T::ColsAtCompileTime == 1;
    };
    */

    template <typename T>
    struct IsVector
    {
        static bool constexpr value = T::RowsAtCompileTime == 1 || T::ColsAtCompileTime == 1;
    };
}