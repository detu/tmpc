#pragma once

#include "Types.hpp"

namespace tmpc :: eigen_adaptor {

template <typename MT>
struct Columns
{
    static size_t constexpr value = MT::ColsAtCompileTime;
};

}