#pragma once

#include "AlignmentFlag.hpp"
#include "PaddingFlag.hpp"
#include "TransposeFlag.hpp"

#include <Eigen/Dense>

namespace tmpc {

template <typename Type, bool AF, bool PF, bool TF = defaultTransposeFlag>
struct CustomVector
:   Eigen::Map<Eigen::Matrix<Type, TF == columnVector ? Eigen::Dynamic : 1, TF == rowVector ? Eigen::Dynamic : 1>>
{
    typedef Eigen::Map<Eigen::Matrix<Type, TF == columnVector ? Eigen::Dynamic : 1, TF == rowVector ? Eigen::Dynamic : 1>> Base;

    CustomVector(Type * data, size_t n)
    :   Base(data, n)
    {        
    }
};

}