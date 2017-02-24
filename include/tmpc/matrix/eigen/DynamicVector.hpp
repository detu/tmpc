#pragma once 

#include "TransposeFlag.hpp"
#include "EigenType.hpp"

#include <Eigen/Dense>

namespace tmpc 
{

    template <typename Type, bool TF>
    using DynamicVectorBase =
        Eigen::Matrix<Type, TF == columnVector ? Eigen::Dynamic : 1, TF == rowVector ? Eigen::Dynamic : 1>;


    template <typename Type, bool TF = defaultTransposeFlag>
    struct DynamicVector
    :   DynamicVectorBase<Type, TF>
    {
        typedef DynamicVectorBase<Type, TF> Base;
        typedef Type ElementType;

        DynamicVector(size_t M)
        :   Base(M)
        {        
        }

        DynamicVector(size_t M, Type const& val)
        :   Base(M)
        {        
            this->setConstant(val);
        }

        template <typename T>
        DynamicVector& operator=(Eigen::MatrixBase<T> const& rhs)
        {
            Base::operator=(rhs);
            return *this;
        }
    };

    template <typename Type, bool TF>
    struct EigenType<DynamicVector<Type, TF>>
    {
        typedef typename DynamicVector<Type, TF>::Base type;
    };

}