#pragma once

#include "Types.hpp"
#include "TransposeFlag.hpp"
#include "VectorAssign.hpp"
#include "EigenType.hpp"

#include <Eigen/Dense>

namespace tmpc {

    template <typename Type, size_t N, bool TF = defaultTransposeFlag>
    using StaticVectorBase = Eigen::Matrix<Type, TF == columnVector ? N : 1, TF == rowVector ? N : 1>;

    template <typename Type, size_t N, bool TF = defaultTransposeFlag>
    struct StaticVector
    :   StaticVectorBase<Type, N, TF>
    ,   Vector<StaticVector<Type, N, TF>, TF>
    {
        typedef StaticVectorBase<Type, N, TF> Base;
        typedef Type ElementType;

        StaticVector()
        {        
        }

        StaticVector(initializer_list<Type> list)
        {
            assign(*this, list);
        }

        template <typename T>
        StaticVector(Eigen::MatrixBase<T> const& rhs)
        :   Base(rhs)
        {        
        }

        StaticVector(Type const& rhs)
        {        
            this->setConstant(rhs);
        }
    };

    template <typename Type, size_t N, bool TF>
    struct EigenType<StaticVector<Type, N, TF>>
    {
        typedef typename StaticVector<Type, N, TF>::Base type;
    };

}