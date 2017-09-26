#pragma once

#include "IsVector.hpp"

namespace tmpc :: eigen_adaptor
{
    template <typename VT>
    struct IsRowVector
    {
        static bool constexpr value = VT::ColsAtCompileTime == 1;
    };
}