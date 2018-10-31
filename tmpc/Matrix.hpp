#pragma once

#include "SizeT.hpp"
#include "matrix/StorageOrder.hpp"
#include "matrix/TransposeFlag.hpp"
#include "matrix/AlignmentFlag.hpp"
#include "matrix/PaddingFlag.hpp"
#include "matrix/BlockMatrixView.hpp"

#include <blaze/Math.h>


namespace tmpc
{
    template <typename Kernel, size_t N, TransposeFlag TF = defaultTransposeFlag>
    using StaticVector = typename Kernel::template StaticVector<N, TF>;

    template <typename Kernel, TransposeFlag TF = defaultTransposeFlag>
    using DynamicVector = typename Kernel::template DynamicVector<TF>;

    template <typename Kernel, size_t M, size_t N, StorageOrder SO = defaultStorageOrder>
    using StaticMatrix = typename Kernel::template StaticMatrix<M, N, SO>;

    template <typename Kernel, StorageOrder SO = defaultStorageOrder>
    using DynamicMatrix = typename Kernel::template DynamicMatrix<SO>;

    template <typename Kernel, AlignmentFlag AF, PaddingFlag PF, StorageOrder SO = defaultStorageOrder>
    using CustomMatrix = typename Kernel::template CustomMatrix<AF, PF, SO>;

    template <typename Kernel, AlignmentFlag AF, PaddingFlag PF, TransposeFlag TF = defaultTransposeFlag>
    using CustomVector = typename Kernel::template CustomVector<AF, PF, TF>;

    template <typename Kernel>
    using IdentityMatrix = typename Kernel::IdentityMatrix;

    template <typename Kernel, typename MT, AlignmentFlag AF = unaligned>
    using Submatrix = typename Kernel::template Submatrix<MT, AF>;

    template <typename Kernel, typename VT, AlignmentFlag AF = unaligned, TransposeFlag TF = defaultTransposeFlag>
    using Subvector = typename Kernel::template Subvector<VT, AF, TF>;

    template <typename Kernel, typename T>
    using Rand = typename Kernel::template Rand<T>;
}