#pragma once

#include <tmpc/matrix/StorageOrder.hpp>
#include "EigenType.hpp"
#include "InitializerList.hpp"
#include "MatrixAssign.hpp"
#include "EigenBase.hpp"
#include "IsEigenMatrix.hpp"

namespace tmpc :: eigen_adaptor 
{
    template <typename MT, StorageOrder SO>
    struct Matrix
    {
        MT& operator~()
        {
            return static_cast<MT&>(*this);
        }

        MT const& operator~() const
        {
            return static_cast<MT const&>(*this);
        }
    };
}

namespace Eigen
{
    template <typename MT>
    inline size_t rows(MatrixBase<MT> const& m)
    {
        return m.rows();
    }

    
    template <typename MT>
    inline size_t columns(MatrixBase<MT> const& m)
    {
        return m.cols();
    }


    template <typename MT>
    inline size_t spacing(MatrixBase<MT> const& m)
    {
        return m.outerStride();
    }


    template <typename MT>
    inline std::pair<size_t, size_t> dimensions(MatrixBase<MT> const& m)
    {
        return std::pair(m.rows(), m.columns());
    }
}