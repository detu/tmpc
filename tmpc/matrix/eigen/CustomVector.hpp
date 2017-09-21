#pragma once

#include <tmpc/matrix/PaddingFlag.hpp>
#include <tmpc/matrix/TransposeFlag.hpp>
#include <tmpc/matrix/AlignmentFlag.hpp>

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor {

template <typename Type, AlignmentFlag AF, PaddingFlag PF, TransposeFlag TF>
using CustomVectorBase =
    Eigen::Map<Eigen::Matrix<Type, TF == columnVector ? Eigen::Dynamic : 1, TF == rowVector ? Eigen::Dynamic : 1>>;

template <typename Type, AlignmentFlag AF, PaddingFlag PF, TransposeFlag TF = defaultTransposeFlag>
struct CustomVector
:   CustomVectorBase<Type, AF, PF, TF>
{
    typedef CustomVectorBase<Type, AF, PF, TF> Base;
    typedef Type ElementType;

    CustomVector(Type * data, size_t n)
    :   Base(data, n)
    {        
    }

    CustomVector& operator=(ElementType const& rhs)
    {
        this->setConstant(rhs);
        return *this;
    }
};

}