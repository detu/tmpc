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
        static_assert(SO1 == blaze::columnMajor);
        static_assert(SO2 == blaze::columnMajor);
        static_assert(SO3 == blaze::columnMajor);

        BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT3 );        
        BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( MT3 );

        if( !isSquare( ~C ) )
            TMPC_THROW_EXCEPTION(std::invalid_argument( "Invalid non-square matrix provided" ));

        size_t const n = rows(C);

        if (columns(A) != n)
            TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid number of columns in matric A"));

        size_t const m = rows(A);

        if( ( !blaze::IsResizable_v<MT2> && ( (~D).rows() != n || (~D).columns() != n ) ) )
            TMPC_THROW_EXCEPTION(std::invalid_argument("Dimensions of fixed size matrix do not match"));

        using Scalar = blaze::ElementType_t<MT3>;

        decltype(auto) l( derestrict( ~D ) );
        resize(l, n, n);

        for (size_t k = 0; k < n; ++k)
        {
            size_t const rs = n - k - 1; // remaining size
        
            auto D01 = submatrix(l, 0, k, k, 1);
            auto D21 = submatrix(l, k + 1, k, rs, 1);
            auto const D10 = submatrix(l, k, 0, 1, k);
            auto const D20 = submatrix(l, k + 1, 0, rs, k);

            // reset(D01);
            
            auto in_C = begin(C, k) + k;
            auto out_D = begin(l, k) + k;

            Scalar x = *in_C++ + dot(column(A, k), column(A, k));

            if (k > 0) 
                x -= sqrNorm(D10);

            if (x <= 0)
                throw std::runtime_error("Unable to continue Cholesky decomposition of a matrix");

            *out_D++ = x = sqrt(x);

            for (size_t i = k + 1; i < n; ++i)
                *out_D++ = *in_C++ + dot(column(A, i), column(A, k));

            for (size_t j = 0; j < k; ++j)
                column(D21, 0) -= (~D)(k, j) * column(D20, j);

            // D21 /= x;
            //
            // NOTE: using scalar/scalar division instead of matrix/scalar division (or multiplication) here
            // because of this issue: https://bitbucket.org/blaze-lib/blaze/issues/288/performance-issue-dmatscalarmultexpr-ctor
            //             
            Scalar const x_inv = 1. / x;
            out_D = begin(l, k) + k + 1;
            for (size_t i = k + 1; i < n; ++i)
                *out_D++ *= x_inv;
        }
    }
}