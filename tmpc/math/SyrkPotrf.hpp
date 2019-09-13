#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/Exception.hpp>

#include <blaze/Math.h>

#include <stdexcept>


namespace tmpc
{
    /// @brief In-place Cholesky decomposition algorithm combined with symmetric rank update.
    ///
    /// Computes Cholesky decomposition D = chol(A^T * A + C).
    ///
    template <typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3>
    inline void syrkPotrf(blaze::DenseMatrix<MT1, SO1> const& A, blaze::DenseMatrix<MT2, SO2> const& C, blaze::DenseMatrix<MT3, SO3>& D)
    {
        BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT3 );        
        BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( MT3 );

        if( !isSquare( ~C ) )
            TMPC_THROW_EXCEPTION(std::invalid_argument( "Invalid non-square matrix provided" ));

        const size_t n = rows(C);

        if (columns(A) != n)
            TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid number of columns in matric A"));

        if( ( !blaze::IsResizable_v<MT2> && ( (~D).rows() != n || (~D).columns() != n ) ) )
            TMPC_THROW_EXCEPTION(std::invalid_argument("Dimensions of fixed size matrix do not match"));

        using Scalar = blaze::ElementType_t<MT3>;

        decltype(auto) l( derestrict( ~D ) );
        resize(l, n, n);

        for (size_t k = 0; k < n; ++k)
        {
            size_t const rs = n - k - 1; // remaining size
        
            auto D01 = submatrix(l, 0, k, k, 1);
            auto D21 = submatrix(l, 0, k, n, 1);
            auto const D10 = submatrix(l, k, 0, 1, k);
            auto const D20 = submatrix(l, 0, 0, n, k);

            reset(D01);
        
            Scalar x = (~C)(k, k) + dot(column(A, k), column(A, k));
            if (k > 0) 
                x -= sqrNorm(D10);

            if (x <= 0)
                throw std::runtime_error("Unable to continue Cholesky decomposition of a matrix");

            l(k, k) = x = sqrt(x);

            // NOTE:
            // Not using matrix-vector multiplication here to prevent Blaze
            // from creating a temporary due to possible aliasing.
            // See this for more details:
            // https://bitbucket.org/blaze-lib/blaze/issues/287/is-there-a-way-to-specify-that-the
            // 
            for (size_t i = k + 1; i < n; ++i)
            {
                D21(i, 0) = (~C)(i, k) + dot(column(A, i), column(A, k)) - row(D20, i) * ctrans(row(D10, 0));

                // NOTE: using scalar/scalar division instead of matrix/scalar division (or multiplication) here
                // because of this issue: https://bitbucket.org/blaze-lib/blaze/issues/288/performance-issue-dmatscalarmultexpr-ctor
                D21(i, 0) /= x;
            }
        }
    }
}