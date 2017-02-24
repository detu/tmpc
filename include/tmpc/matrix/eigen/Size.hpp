#pragma once

#include "Types.hpp"

#include <Eigen/Dense>

#include <type_traits>

namespace tmpc
{
    template <typename T, typename Enable = void>
    struct Size;

    template <typename T>
    struct Size<T, 
        std::enable_if_t<IsVector<T>::value>
    >
    {
        static size_t constexpr value = T::SizeAtCompileTime;
    };
}