/// \brief Some useful math functions.

#pragma once

#include <blaze/Math.h>

#include <limits>


namespace tmpc
{
    template <typename T>
    inline T constexpr inf()
    {
        return std::numeric_limits<T>::infinity();
    }


    template <typename T>
    inline T constexpr sNaN()
    {
        return std::numeric_limits<T>::signaling_NaN();
    }
}


namespace blaze
{
    /// @brief No-resize adaptor for vector assignment
    ///
    /// See https://bitbucket.org/blaze-lib/blaze/issues/93/matrix-vector-assignment-without-resizing
    template <typename VT, bool TF>
    inline auto noresize(Vector<VT, TF>& v)
    {
        return subvector(v, 0, size(v));
    }


    /// @brief No-resize adaptor for matrix assignment
    ///
    /// See https://bitbucket.org/blaze-lib/blaze/issues/93/matrix-vector-assignment-without-resizing
    template <typename MT, bool SO>
    inline auto noresize(Matrix<MT, SO>& m)
    {
        return submatrix(m, 0, 0, rows(m), columns(m));
    }
}