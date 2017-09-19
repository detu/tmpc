#pragma once

#include <tmpc/matrix/TransposeFlag.hpp>
#include <tmpc/matrix/StorageOrder.hpp>
#include "Types.hpp"

namespace tmpc :: eigen_adaptor
{
    template <typename Type, size_t N, TransposeFlag TF>
    struct StaticVector;

    template <typename Type, TransposeFlag TF>
    struct DynamicVector;

    template <typename VT, AlignmentFlag AF, TransposeFlag TF>
    struct Subvector;

    template <typename Type, size_t M, size_t N, StorageOrder SO>
    struct StaticMatrix;

    template <typename Type, StorageOrder SO>
    struct DynamicMatrix;

    template <typename MT, AlignmentFlag AF>
    struct Submatrix;
}