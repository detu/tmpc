#pragma once

#include <tmpc/Math.hpp>
#include <tmpc/Exception.hpp>


namespace tmpc
{
    /// @brief Solve the equation A*x=b with a triangular matrix A
    ///
    template <
        typename MT, // Type of the system matrix
        bool SO, // Storage order of the system matrix
        typename VT1, // Type of the right-hand side vector
        typename VT2 // Type of the result vector
    >
    inline void trsv(blaze::DenseMatrix<MT, SO> const& A, 
        blaze::DenseVector<VT1, blaze::columnVector> const& b, 
        blaze::DenseVector<VT2, blaze::columnVector>& x)
    {
        BLAZE_CONSTRAINT_MUST_BE_TRIANGULAR_MATRIX_TYPE( MT );
        BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

        BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
        BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
        BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT1 );

        BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );
        BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT2 );
        BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT2 );

        if (!isSquare(~A))
            TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid non-square matrix provided"));

        auto const n = rows(A);

        if (size(b) != n)
            TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid right-hand side vector provided"));
            
        if (size(x) != n)
            TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid result vector provided"));
        

        if (blaze::IsLower_v<MT>)
        {
            // Lower matrix algorithm

            if (blaze::IsRowMajorMatrix_v<MT>)
            {
                // Row-major algorithm
                for (size_t i = 0; i < n; ++i)
                    (~x)[i] = ((~b)[i] - subvector(row(~A, i), 0, i) * subvector(~x, 0, i)) / (~A)(i, i);
            }
            else
            {
                // Column-major algorithm
                ~x = ~b;

                for (size_t i = 0; i < n; ++i)
                {
                    auto const x_i = (~x)[i] /= (~A)(i, i);
                    subvector(~x, i + 1, n - i - 1) -= subvector(column(~A, i), i + 1, n - i - 1) * x_i;
                }
            }
        }
        else
        {
            // Upper matrix algorithm

            if (blaze::IsRowMajorMatrix_v<MT>)
            {
                // Row-major algorithm
                for (size_t i = n; i-- > 0; )
                    (~x)[i] = ((~b)[i] - subvector(row(~A, i), i + 1, n - i - 1) * subvector(~x, i + 1, n - i - 1)) / (~A)(i, i);
            }
            else
            {
                // Column-major algorithm
                ~x = ~b;

                for (size_t i = n; i-- > 0; )
                {
                    auto const x_i = (~x)[i] /= (~A)(i, i);
                    subvector(~x, 0, i) -= subvector(column(~A, i), 0, i) * x_i;
                }
            }
        }
    }


    /// @brief Solve the equation x*A=b with a triangular matrix A and a row vector x.
    ///
    template <
        typename VT1, // Type of the right-hand side vector
        typename MT, // Type of the system matrix
        bool SO, // Storage order of the system matrix
        typename VT2 // Type of the result vector
    >
    inline void trsv(
        blaze::DenseVector<VT1, blaze::rowVector> const& b,
        blaze::DenseMatrix<MT, SO> const& A,
        blaze::DenseVector<VT2, blaze::rowVector>& x)
    {
        BLAZE_CONSTRAINT_MUST_BE_TRIANGULAR_MATRIX_TYPE( MT );
        BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

        BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
        BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
        BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT1 );

        BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );
        BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT2 );
        BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT2 );

        if (!isSquare(~A))
            TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid non-square matrix provided"));

        auto const n = rows(A);

        if (size(b) != n)
            TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid right-hand side vector provided"));
            
        if (size(x) != n)
            TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid result vector provided"));
        

        if (blaze::IsUpper_v<MT>)
        {
            // Upper matrix algorithm

            if (blaze::IsColumnMajorMatrix_v<MT>)
            {
                // Column-major algorithm
                for (size_t i = 0; i < n; ++i)
                    (~x)[i] = ((~b)[i] - subvector(~x, 0, i) * subvector(column(~A, i), 0, i)) / (~A)(i, i);
            }
            else
            {
                // Row-major algorithm
                ~x = ~b;

                for (size_t i = 0; i < n; ++i)
                {
                    auto const x_i = (~x)[i] /= (~A)(i, i);
                    subvector(~x, i + 1, n - i - 1) -= x_i * subvector(row(~A, i), i + 1, n - i - 1);
                }
            }
        }
        else
        {
            // Lower matrix algorithm

            if (blaze::IsColumnMajorMatrix_v<MT>)
            {
                // Column-major algorithm
                for (size_t i = n; i-- > 0; )
                    (~x)[i] = ((~b)[i] - subvector(~x, i + 1, n - i - 1) * subvector(column(~A, i), i + 1, n - i - 1)) / (~A)(i, i);
            }
            else
            {
                // Row-major algorithm
                ~x = ~b;

                for (size_t i = n; i-- > 0; )
                {
                    auto const x_i = (~x)[i] /= (~A)(i, i);
                    subvector(~x, 0, i) -= x_i * subvector(row(~A, i), 0, i);
                }
            }
        }
    }
}
