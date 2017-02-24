#pragma once

#include "TransposeFlag.hpp"
#include "StorageOrder.hpp"
#include "Types.hpp"

namespace tmpc
{
    template <typename Type, bool SO>
    struct DynamicMatrix;

    template <typename Type, size_t M, size_t N, bool SO>
    struct StaticMatrix;

    template <typename MT, bool AF>
    struct Submatrix;
}