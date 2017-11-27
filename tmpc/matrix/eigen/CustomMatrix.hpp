#pragma once

#include <tmpc/matrix/PaddingFlag.hpp>
#include <tmpc/matrix/StorageOrder.hpp>
#include <tmpc/matrix/AlignmentFlag.hpp>

#include "Eigen.hpp"
#include "Types.hpp"
#include "Matrix.hpp"

#include <stdexcept>

namespace tmpc :: eigen_adaptor 
{
    template <typename Type, StorageOrder SO>
    using CustomMatrixBaseUnalignedUnpadded =
        Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor>>;

    template <typename Type, AlignmentFlag AF, PaddingFlag PF, StorageOrder SO = defaultStorageOrder>
    struct CustomMatrix;

    template <typename Type, StorageOrder SO>
    struct CustomMatrix<Type, unaligned, unpadded, SO> 
    :   Matrix<CustomMatrix<Type, unaligned, unpadded, SO>, SO>
    ,   CustomMatrixBaseUnalignedUnpadded<Type, SO>
    {
        using EigenBase = CustomMatrixBaseUnalignedUnpadded<Type, SO>;
        using ElementType = Type;

        CustomMatrix(ElementType * ptr, size_t m, size_t n)
        :   EigenBase {ptr, static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(n)}
        {        
        }

        CustomMatrix& operator=(ElementType const& rhs)
        {
            this->setConstant(rhs);
            return *this;
        }

        void reset(ElementType * ptr, size_t m, size_t n)
        {
            if (this->rows() != m || this->cols() != n)
                throw std::runtime_error("Impossible to change the size of CustomMatrix<BlazeKernel> in reset()");

            this->m_data = ptr;        
        }
    };
}