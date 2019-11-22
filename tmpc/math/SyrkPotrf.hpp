#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/Exception.hpp>

#include <blazefeo/math/DynamicPanelMatrix.hpp>
#include <blazefeo/math/panel/Potrf.hpp>

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

        l = declsym(~C + trans(~A) * ~A);

        static_assert(SO3 == blaze::columnMajor);
        blazefeo::DynamicPanelMatrix<Scalar, SO3> D1(n, n);
        D1 = l;

        potrf(D1, D1);
        D1.unpackLower(data(l), spacing(l));

        for (size_t j = 0; j < n; ++j)
            reset(subvector(column(l, j), 0, j));
    }
}