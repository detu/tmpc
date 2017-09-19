#pragma once

#include <tmpc/matrix/StorageOrder.hpp>
#include "EigenType.hpp"
#include "InitializerList.hpp"
#include "MatrixAssign.hpp"
#include "EigenBase.hpp"
#include "IsEigenMatrix.hpp"

namespace tmpc :: eigen_adaptor 
{
    struct MatrixTag {};

    template <typename MT, StorageOrder SO>
    struct Matrix : MatrixTag
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

    template <typename MT>
    std::enable_if_t<IsEigenMatrix<MT>::value, size_t> rows(MT const& m)
    {
        return m.rows();
    }

    template <typename MT, StorageOrder SO>
    size_t rows(Matrix<MT, SO> const& m)
    {
        return (~m).rows();
    }

    template <typename MT>
    std::enable_if_t<IsEigenMatrix<MT>::value, size_t> columns(MT const& m)
    {
        return m.cols();
    }

    template <typename MT, StorageOrder SO>
    size_t columns(Matrix<MT, SO> const& m)
    {
        return (~m).cols();
    }
}