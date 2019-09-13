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
            // auto D01 = submatrix(l, 0, k, k, 1);
            // reset(D01);
        
            auto w_col = column(l, k);
            auto const w_mat = submatrix(l, 0, 0, n, k);
            auto const w_row = row(w_mat, k);            
        
            Scalar x = (~C)(k, k) + dot(column(A, k), column(A, k)) - sqrNorm(w_row);

            if (x <= 0)
                TMPC_THROW_EXCEPTION(std::runtime_error("Unable to continue Cholesky decomposition of a matrix"));

            l(k, k) = x = sqrt(x);

            // NOTE:
            // Not using matrix-vector multiplication here to prevent Blaze
            // from creating a temporary due to possible aliasing.
            // See this for more details:
            // https://bitbucket.org/blaze-lib/blaze/issues/287/is-there-a-way-to-specify-that-the
            // 
            for (size_t i = k + 1; i < n; ++i)
                w_col[i] = ((~C)(i, k) + dot(column(A, i), column(A, k)) - dot(row(w_mat, i), w_row)) / x;
        }
    }
}