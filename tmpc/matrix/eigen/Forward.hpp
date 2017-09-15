#pragma once

#include "TransposeFlag.hpp"
#include "StorageOrder.hpp"
#include "Types.hpp"

namespace tmpc
{
    template <typename Type, size_t N, bool TF>
    struct StaticVector;

    template <typename Type, bool TF>
    struct DynamicVector;

    template <typename VT, bool AF, bool TF>
    struct Subvector;

    template <typename Type, size_t M, size_t N, bool SO>
    struct StaticMatrix;

    template <typename Type, bool SO>
    struct DynamicMatrix;

    template <typename MT, bool AF>
    struct Submatrix;
}