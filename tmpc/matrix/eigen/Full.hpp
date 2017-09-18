#pragma once

#include "Subvector.hpp"
#include "Submatrix.hpp"
#include "IsVector.hpp"

#include <type_traits>

namespace tmpc :: eigen_adaptor
{
    template <typename VT>
    std::enable_if_t<
        IsVector<VT>::value, 
        Subvector<VT, unaligned>
    > full(VT& v)
    {
        return subvector(v, 0, size(v));
    }

    template <typename MT>
    Submatrix<MT, unaligned> full(Eigen::MatrixBase<MT>& m)
    {
        return submatrix(m, 0, 0, m.rows(), m.cols());
    }
}