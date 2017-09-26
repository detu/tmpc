#pragma once

#include "Types.hpp"
#include "IsVector.hpp"

#include "Eigen.hpp"

#include <type_traits>

namespace tmpc :: eigen_adaptor
{
    template <typename T, typename Enable = void>
    struct Size;

    template <typename T>
    struct Size<T, 
        std::enable_if_t<IsVector<T>::value>
    >
    {
        static auto constexpr value = T::SizeAtCompileTime;
    };
}