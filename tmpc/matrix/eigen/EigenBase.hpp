#pragma once

//#include "Forward.hpp"
#include "IsEigenMatrix.hpp"

#include "Eigen.hpp"

namespace tmpc
{
    template <typename T>
    struct EigenBaseSelector;

    template <typename T>
    using EigenBase = typename EigenBaseSelector<T>::type;

    template <typename T>
    struct EigenBaseSelector //<T, std::enable_if_t<IsEigenMatrix<T>::value>>
    {
        using type = T;
    };

    template <typename T>
    struct EigenBaseSelector<T const>
    {
        using type = EigenBase<T> const;
    };
}