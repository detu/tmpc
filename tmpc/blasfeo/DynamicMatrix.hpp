#pragma once

#include <tmpc/blasfeo/BlasfeoMatrix.hpp>
#include <tmpc/blasfeo/Matrix.hpp>
#include <tmpc/blasfeo/Memory.hpp>


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


        DynamicMatrix(size_t m, size_t n)
        :   data_(detail::alignedAlloc<Real>(memsize_mat<Real>(m, n)))
        {
            create_mat(m, n, *this, data_.get());
        }


        /// @brief Resize the matrix to new size and re-allocate memory if needed.
        ///
        /// Does not preserve matrix elements if reallocation occurs.
        void resize(size_t m, size_t n)
        {
            auto const required_capacity = memsize_mat<Real>(m, n);

            if (required_capacity > capacity_)
            {
                data_.reset(detail::alignedAlloc<Real>(required_capacity));
                capacity_ = required_capacity;
            }

            create_mat(m, n, *this, data_.get());
        }


    private:
        std::unique_ptr<Real[], detail::AlignedDeleter> data_;

        /// @brief Actually allocated bytes pointed by data_
        size_t capacity_ = 0;
    };
}