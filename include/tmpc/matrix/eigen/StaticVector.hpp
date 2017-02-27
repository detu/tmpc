#pragma once

#include "Types.hpp"
#include "TransposeFlag.hpp"
#include "VectorAssign.hpp"
//#include "EigenType.hpp"
#include "EigenBase.hpp"

namespace tmpc 
{
	template <typename Type, size_t N, bool TF>
    struct StaticVector;
	
    template <typename Type, size_t N, bool TF>
    struct EigenBaseSelector<StaticVector<Type, N, TF>>
    {
        using type = Eigen::Matrix<Type, TF == columnVector ? N : 1, TF == rowVector ? N : 1>;
    };

    template <typename Type, size_t N, bool TF = defaultTransposeFlag>
    struct StaticVector
    :   EigenBase<StaticVector<Type, N, TF>>
    ,   Vector<StaticVector<Type, N, TF>, TF>
    {
        typedef EigenBase<StaticVector<Type, N, TF>> Base;
        typedef Type ElementType;

        StaticVector()
        {        
        }

        explicit StaticVector(initializer_list<Type> list)
        {
            assign(*this, list);
        }

        template <typename T>
        StaticVector(Eigen::MatrixBase<T> const& rhs)
        :   Base(rhs)
        {        
        }

        explicit StaticVector(ElementType const& rhs)
        {        
            this->setConstant(rhs);
        }

        StaticVector& operator=(ElementType const& rhs)
        {
            this->setConstant(rhs);
            return *this;
        }

        template <typename T>
        StaticVector& operator=(Eigen::MatrixBase<T> const& rhs)
        {
            Base::operator=(rhs);
            return *this;
        }
    };

    /*
    template <typename Type, size_t N, bool TF>
    struct EigenType<StaticVector<Type, N, TF>>
    {
        typedef typename StaticVector<Type, N, TF>::Base type;
    };
    */
}