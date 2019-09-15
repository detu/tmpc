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
            if constexpr(!blaze::IsLower_v<MT3>)
                reset(submatrix(l, 0, k, k, 1));
            
            auto in_C = begin(C, k) + k;
            auto out_D = begin(l, k) + k;
            
        #if 1
            for (size_t j = k; j < n; ++j)
                *out_D++ = *in_C++ + dot(column(A, j), column(A, k));
        #else
            for (size_t j = k; j < n; ++j, ++in_C, ++out_D)
            {
                auto in_Aj = begin(~A, j);
                auto in_Ak = begin(~A, k);

                SIMDType simd_s = blaze::set(0.);
                size_t jj = 0;

                for ( ; jj + SIMDSIZE <= m; jj += SIMDSIZE)
                {
                    simd_s += in_Aj.load() * in_Ak.load();
                    in_Aj += SIMDSIZE;
                    in_Ak += SIMDSIZE;
                }

                Scalar s = sum(simd_s);

                for ( ; jj < m; ++jj)
                    s += *in_Aj++ * *in_Ak++;

                *out_D = *in_C + s;
            }
        #endif


        #if 1
            // TODO: rewrite this loop as matrix-vector multiplication when this issue is fixed:
            // https://bitbucket.org/blaze-lib/blaze/issues/287/is-there-a-way-to-specify-that-the
            for (size_t j = 0; j < k; ++j)
                // Blaze detects that the columns are different and therefore there is no aliasing.
                column(D21, 0) -= (~D)(k, j) * column(D20, j);
        
        #else
            for (size_t j = 0; j < k; ++j)
            {
                auto in = begin(l, j) + k;
                auto out = begin(l, k) + k;

                Scalar const D_kj = *in;
                SIMDType const simd_D_kj = blaze::set(D_kj);
                size_t i = k;

                bool constexpr remainder = true;
                size_t const ipos( ( remainder )?( n & size_t(-SIMDSIZE) ):( n ) );

                for ( ; i < ipos; i += SIMDSIZE)
                {
                    out.storeu(out.loadu() - in.loadu() * simd_D_kj);
                    in += SIMDSIZE;
                    out += SIMDSIZE;
                }

                for ( ; i < n; ++i)
                    *out++ -= *in++ * D_kj;
            }
        #endif

            out_D = begin(l, k) + k;
            Scalar x = *out_D;
            if (x <= 0)
                throw std::runtime_error("Unable to continue Cholesky decomposition of a matrix");

            *out_D++ = x = sqrt(x);

            // D21 /= x;
            //
            // NOTE 1: using scalar/scalar division instead of matrix/scalar division (or multiplication) here
            // because of this issue: https://bitbucket.org/blaze-lib/blaze/issues/288/performance-issue-dmatscalarmultexpr-ctor
            //
            // NOTE 2: from looking from assembly, this loop is not vectorized.
            //
            Scalar const x_inv = 1. / x;
            for (size_t i = k + 1; i < n; ++i)
                *out_D++ *= x_inv;
        }
    }
}