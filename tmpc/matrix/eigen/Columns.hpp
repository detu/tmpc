#pragma once

#include "Types.hpp"

namespace tmpc {

template <typename MT>
struct Columns
{
    static size_t constexpr value = MT::ColsAtCompileTime;
};

}