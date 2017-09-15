#pragma once

#include "Types.hpp"

namespace tmpc {

template <typename MT>
struct Rows
{
    static size_t constexpr value = MT::RowsAtCompileTime;
};

}