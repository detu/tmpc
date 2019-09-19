#pragma once

#include <tmpc/blasfeo/BlasfeoMatrix.hpp>
#include <tmpc/blasfeo/Matrix.hpp>
#include <tmpc/blasfeo/Memory.hpp>

#include <blaze/Math.h>


namespace tmpc :: blasfeo
{
    /// @brief BLASFEO matrix managing its own memory pool
    ///
    template <typename Real>
    class DynamicMatrix
    :   public BlasfeoMatrix_t<Real>
    ,   public Matrix<DynamicMatrix<Real>>
    {
    public:
        using ElementType = Real;

    
        /// \brief Create a 0-by-0 matrix.
        DynamicMatrix()
        {
            create_mat(0, 0, *this, nullptr);
        }


        /// \brief Create a matrix of given size.
        DynamicMatrix(size_t m, size_t n)
        :   data_(memsize_mat<Real>(m, n))
        {
            create_mat(m, n, *this, data_.get());
        }


        /// \brief Create a copy of a Blaze dense column-major matrix.
        template <typename MT>
        DynamicMatrix(blaze::DenseMatrix<MT, blaze::columnMajor> const& rhs)
        :   DynamicMatrix(rows(rhs), columns(rhs))
        {
            *this = rhs;
        }


        /// @brief Resize the matrix to new size and re-allocate memory if needed.
        ///
        /// Does not preserve matrix elements if reallocation occurs.
        void resize(size_t m, size_t n)
        {
            if (m != rows(*this) || n != columns(*this))
            {
                auto const required_capacity = memsize_mat<Real>(m, n);

                if (required_capacity > capacity_)
                {
                    data_ = detail::AlignedStorage(required_capacity);
                    capacity_ = required_capacity;
                }

                create_mat(m, n, *this, data_.get());
            }
        }


        /// @brief Assign Blaze column-major dense matrix to BLASFEO matrix.
        ///
        /// The BLASFEO matrix is resized if needed.
        template <typename MT>
        DynamicMatrix& operator=(blaze::DenseMatrix<MT, blaze::columnMajor> const& rhs)
        {
            auto const m = rows(rhs);
            auto const n = columns(rhs);

            resize(m, n);
            pack_mat(m, n, data(rhs), spacing(rhs), *this, 0, 0);

            return *this;
        }


    private:
        detail::AlignedStorage data_;

        /// @brief Actually allocated bytes pointed by data_
        size_t capacity_ = 0;
    };
}