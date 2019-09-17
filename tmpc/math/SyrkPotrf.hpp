#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/Exception.hpp>

#include <blaze/Math.h>


namespace tmpc
{
    /// @brief In-place Cholesky decomposition algorithm combined with symmetric rank update.
    ///
    /// Computes Cholesky decomposition D = chol(A^T * A + C).
    ///
    template <typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3>
    inline void syrkPotrf(blaze::DenseMatrix<MT1, SO1> const& A, blaze::DenseMatrix<MT2, SO2> const& C, blaze::DenseMatrix<MT3, SO3>& D)
    {
        BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE(MT1);
        BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE(MT2);

        BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT3 );        
        BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( MT3 );
        BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE(MT3);
        
        if( !isSquare( ~C ) )
            TMPC_THROW_EXCEPTION(std::invalid_argument( "Invalid non-square matrix provided" ));

        size_t const n = rows(C);

        if (columns(A) != n)
            TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid number of columns in matric A"));

        size_t const m = rows(A);

        if( ( !blaze::IsResizable_v<MT2> && ( (~D).rows() != n || (~D).columns() != n ) ) )
            TMPC_THROW_EXCEPTION(std::invalid_argument("Dimensions of fixed size matrix do not match"));

        using Scalar = blaze::ElementType_t<MT3>;
        BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE(blaze::ElementType_t<MT1>, Scalar);
        BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE(blaze::ElementType_t<MT2>, Scalar);

        size_t constexpr SIMDSIZE = blaze::SIMDTrait<Scalar>::size;
        using SIMDType = blaze::SIMDTrait_t<Scalar>;
        bool constexpr use_simd = false;

        decltype(auto) l( derestrict( ~D ) );
        resize(l, n, n);

        for (size_t k = 0; k < n; ++k)
        {
            size_t const rs = n - k; // remaining size
        
            auto D21 = submatrix(l, k, k, rs, 1);
            auto const D20 = submatrix(l, k, 0, rs, k);

            // Set upper-triangular elements to 0 if the result matrix is not restricted lower
            if constexpr (!blaze::IsLower_v<MT3>)
                reset(submatrix(l, 0, k, k, 1));
            
        #if 1
            for (size_t j = k; j < n; ++j)
                l(j, k) = (~C)(j, k) + dot(column(A, j), column(A, k));
        #else
            for (size_t j = k; j < n; ++j)
            {
                SIMDType xmm0, xmm1, xmm2, xmm3;
                size_t jj = 0;

                for ( ; jj + 4 * SIMDSIZE <= m; jj += 4 * SIMDSIZE)
                {
                    xmm0 = xmm0 + (~A).load(jj, j) * (~A).load(jj, k);
                    xmm1 = xmm1 + (~A).load(jj + 1 * SIMDSIZE, j) * (~A).load(jj + 1 * SIMDSIZE, k);
                    xmm2 = xmm2 + (~A).load(jj + 2 * SIMDSIZE, j) * (~A).load(jj + 2 * SIMDSIZE, k);
                    xmm3 = xmm3 + (~A).load(jj + 3 * SIMDSIZE, j) * (~A).load(jj + 3 * SIMDSIZE, k);
                }

                for ( ; jj + 2 * SIMDSIZE <= m; jj += 2 * SIMDSIZE)
                {
                    xmm0 = xmm0 + (~A).load(jj, j) * (~A).load(jj, k);
                    xmm1 = xmm1 + (~A).load(jj + SIMDSIZE, j) * (~A).load(jj + SIMDSIZE, k);
                }

                for ( ; jj + SIMDSIZE <= m; jj += SIMDSIZE)
                    xmm0 = xmm0 + (~A).load(jj, j) * (~A).load(jj, k);

                Scalar s = sum(xmm0 + xmm1 + xmm2 + xmm3);

                for ( ; jj < m; ++jj)
                    s += (~A)(jj, j) * (~A)(jj, k);

                *out_D = *in_C + s;
            }
        #endif


        #if 1
            // TODO: rewrite this loop as matrix-vector multiplication when this issue is fixed:
            // https://bitbucket.org/blaze-lib/blaze/issues/287/is-there-a-way-to-specify-that-the
            //
            // submatrix(l, k, k, rs, 1) -= submatrix(l, k, 0, rs, k) * trans(submatrix(l, k, 0, 1, k));
            //
            for (size_t j = 0; j < k; ++j)
                // NOTE 1:
                // Blaze detects that the columns are different and therefore there is no aliasing.
                //
                // NOTE 2:
                // "unchecked" on the left-hand side of the assignment seems to make the code faster,
                // whereas on the right side it slows it dows.
                //
                // NOTE: regarding checked vs unchecked performance, see this issue:
                // https://bitbucket.org/blaze-lib/blaze/issues/291/unchecked-column-of-a-checked-submatrix-is
                //
                column(D21, 0) -= (~D)(k, j) * column(D20, j);
        #else
            for (size_t j = 0; j < k; ++j)
            {
                Scalar const D_kj = (~D)(k, j);
                SIMDType const simd_D_kj = blaze::set(D_kj);
                size_t i = k;

                // bool constexpr remainder = true;
                // size_t const ipos( ( remainder )?( n & size_t(-SIMDSIZE) ):( n ) );

                for ( ; i % SIMDSIZE != 0; ++i)
                    l(i, k) -= l(i, j) * D_kj;

                for ( ; i + 4 * SIMDSIZE <= n; i += 4 * SIMDSIZE)
                {
                    l.storea(i, k, l.loada(i, k) - l.loada(i, j) * simd_D_kj);
                    l.storea(i + 1 * SIMDSIZE, k, l.loada(i + 1 * SIMDSIZE, k) - l.loada(i + 1 * SIMDSIZE, j) * simd_D_kj);
                    l.storea(i + 2 * SIMDSIZE, k, l.loada(i + 2 * SIMDSIZE, k) - l.loada(i + 2 * SIMDSIZE, j) * simd_D_kj);
                    l.storea(i + 3 * SIMDSIZE, k, l.loada(i + 3 * SIMDSIZE, k) - l.loada(i + 3 * SIMDSIZE, j) * simd_D_kj);
                }

                for ( ; i + 2 * SIMDSIZE <= n; i += 2 * SIMDSIZE)
                {
                    l.storea(i, k, l.loada(i, k) - l.loada(i, j) * simd_D_kj);
                    l.storea(i + SIMDSIZE, k, l.loada(i + SIMDSIZE, k) - l.loada(i + SIMDSIZE, j) * simd_D_kj);
                }

                for ( ; i + SIMDSIZE <= n; i += SIMDSIZE)
                    l.storea(i, k, l.loada(i, k) - l.loada(i, j) * simd_D_kj);

                for ( ; i < n; ++i)
                    l(i, k) -= l(i, j) * D_kj;
            }
        #endif

            Scalar x = l(k, k);
            if (x <= 0)
                throw std::runtime_error("Unable to continue Cholesky decomposition of a matrix");

            // D21 /= x;
            //
            // NOTE 1: using scalar/scalar division instead of matrix/scalar division (or multiplication) here
            // because of this issue: https://bitbucket.org/blaze-lib/blaze/issues/288/performance-issue-dmatscalarmultexpr-ctor
            //
            // NOTE 2: from looking from assembly, this loop is not vectorized.
            //
            // NOTE 3: don't do l(k, k) = sqrt(x) to avoid transferring a single value to memory.
            // Instead, assign it as a part of vectorized multiply operation.
            //
            Scalar const x_inv = 1. / sqrt(x);
            SIMDType simd_x_inv = blaze::set(x_inv);
            size_t i = k;

            for ( ; i % SIMDSIZE != 0; ++i)
                l(i, k) *= x_inv;

            for ( ; i + 4 * SIMDSIZE <= n; i += 4 * SIMDSIZE)
            {
                l.storea(i, k, l.loada(i, k) * simd_x_inv);
                l.storea(i + 1 * SIMDSIZE, k, l.loada(i + 1 * SIMDSIZE, k) * simd_x_inv);
                l.storea(i + 2 * SIMDSIZE, k, l.loada(i + 2 * SIMDSIZE, k) * simd_x_inv);
                l.storea(i + 3 * SIMDSIZE, k, l.loada(i + 3 * SIMDSIZE, k) * simd_x_inv);
            }

            for ( ; i + 2 * SIMDSIZE <= n; i += 2 * SIMDSIZE)
            {
                l.storea(i, k, l.loada(i, k) * simd_x_inv);
                l.storea(i + SIMDSIZE, k, l.loada(i + SIMDSIZE, k) * simd_x_inv);
            }

            for ( ; i < n; ++i)
                l(i, k) *= x_inv;
        }
    }
}