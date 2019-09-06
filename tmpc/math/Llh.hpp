#pragma once

#include <tmpc/SizeT.hpp>

#include <blaze/Math.h>

#include <stdexcept>


namespace tmpc
{
    /// @brief Cholesky decomposition algorithm.
    ///
    template <typename MT1, bool SO1, typename MT2, bool SO2>
    inline void llh(blaze::DenseMatrix<MT1, SO1> const& A, blaze::DenseMatrix<MT2, SO2>& L)
    {
        BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
        BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( blaze::ElementType_t<MT1> );

        BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT2 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT2 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( MT2 );
        BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( blaze::ElementType_t<MT2> );

        if( !isSquare( ~A ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
        }

        const size_t n( (~A).rows() );

        if( ( !blaze::IsResizable_v<MT2> && ( (~L).rows() != n || (~L).columns() != n ) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Dimensions of fixed size matrix do not match" );
        }

        using Scalar = blaze::ElementType_t<MT1>;

        // Copy A to L
        decltype(auto) l( derestrict( ~L ) );
        l = ~A;

        for (size_t k = 0; k < n; ++k)
        {
            size_t const rs = n - k - 1; // remaining size
        
            auto A01 = submatrix(l, 0, k, k, 1);
            auto A21 = submatrix(l, k + 1, k, rs, 1);
            auto const A10 = submatrix(l, k, 0, 1, k);
            auto const A20 = submatrix(l, k + 1, 0, rs, k);

            reset(A01);
        
            Scalar x = l(k, k);
            if (k > 0) 
                x -= sqrNorm(A10);

            if (x <= 0)
                throw std::runtime_error("Unable to continue Cholesky decomposition of a matrix");

            l(k, k) = x = sqrt(x);

            if (k > 0 && rs > 0) 
                A21 -= A20 * ctrans(A10);

            if (rs > 0) 
                A21 /= x;
        }
    }
}