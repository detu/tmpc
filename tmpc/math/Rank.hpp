#pragma once

#include <blaze/Math.h>

#include <algorithm>
#include <limits>


namespace tmpc
{
    /// @brief Computes rank of a matrix
    ///
    /// TODO: needs to be replaced by a corresponding blaze function once it is implemented:
    /// https://bitbucket.org/blaze-lib/blaze/issues/264
    ///
    /// @param A the whose rank is to be computed.
    /// @return Rank of matrix \a A.
    template <typename MT, bool SO>
    inline blaze::size_t rank(blaze::Matrix<MT, SO> const& A)
    {
        using ET = blaze::ElementType_t<MT>;
        blaze::DynamicVector<ET> s;

        // Compute singular values of A
        svd(~A, s);

        // Count singular values greater than the tolerance.
        // The tolerance is computed as max(size(A))*eps(norm(A)), 
        // in accordance with how the MATLAB rank() function works:
        // https://www.mathworks.com/help/matlab/ref/rank.html
        ET const tol = std::max(rows(A), columns(A)) * std::numeric_limits<ET>::epsilon() * max(s);
        return std::count_if(begin(s), end(s), [tol] (auto val) { return abs(val) > tol; });
    }
}