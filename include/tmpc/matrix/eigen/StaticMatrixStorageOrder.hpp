#pragma once

#include "Types.hpp"
#include "StorageOrder.hpp"

#include <Eigen/Dense>

namespace tmpc
{
    template <size_t M, size_t N, bool SO>
    struct StaticMatrixStorageOrder
    {
        static bool constexpr value = 
            (M > 1 && N == 1) 
            ? Eigen::ColMajor 
            : (SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor);
    };
}