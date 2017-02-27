#pragma once

#include "AlignmentFlag.hpp"
#include "PaddingFlag.hpp"
#include "TransposeFlag.hpp"

#include "Eigen.hpp"

namespace tmpc {

template <typename Type, bool AF, bool PF, bool TF>
using CustomVectorBase =
    Eigen::Map<Eigen::Matrix<Type, TF == columnVector ? Eigen::Dynamic : 1, TF == rowVector ? Eigen::Dynamic : 1>>;

template <typename Type, bool AF, bool PF, bool TF = defaultTransposeFlag>
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