#pragma once

#include <blaze/Math.h>


namespace tmpc
{
    /// @brief Matrix exponential.
    ///
    template <typename MT, bool SO>
    [[deprecated("Replaced by blaze::matexp()")]]
    inline decltype(auto) expm(blaze::DenseMatrix<MT, SO> const& m)
    {
        return blaze::matexp(m);
    }
}