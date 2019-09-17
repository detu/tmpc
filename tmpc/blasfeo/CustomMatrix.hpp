#pragma once

#include <tmpc/blasfeo/BlasfeoMatrix.hpp>
#include <tmpc/blasfeo/Matrix.hpp>


namespace tmpc :: blasfeo
{
    /// @brief BLASFEO matrix using preallocated memory array
    ///
    template <typename Real>
    class CustomMatrix
    :   public BlasfeoMatrix_t<Real>
    ,   public Matrix<CustomMatrix<Real>>
    {
    public:
        using ElementType = Real;


        /// \brief Create a 0-by-0 matrix.
        CustomMatrix()
        {
            create_mat(0, 0, this, nullptr);
        }

        
        CustomMatrix(Real * data, size_t m, size_t n)
        {
            create_mat(m, n, *this, data);
        }


        /// \brief Set new pointer and dimensions.
        /// Use with care!
        void reset(Real * data, size_t m, size_t n)
        {
            create_mat(m, n, *this, data);
        }
    };
}