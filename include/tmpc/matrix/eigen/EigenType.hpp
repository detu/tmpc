#pragma once

#include <Eigen/Dense>

#include <type_traits>

namespace tmpc
{
    template <typename T, typename Enable = void>
    struct EigenType;

    template <typename T>
    struct EigenType<
        T,
        std::enable_if_t< 
            std::is_base_of<Eigen::EigenBase<T>, T>::value
        >
    >
    {
        typedef T type;
    };
}