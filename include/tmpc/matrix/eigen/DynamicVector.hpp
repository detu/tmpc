#pragma once 

#include "TransposeFlag.hpp"
//#include "EigenType.hpp"
#include "Vector.hpp"
#include "EigenBase.hpp"
#include "Rand.hpp"

namespace tmpc 
{
    template <typename Type, bool TF>
    struct DynamicVector;

    template <typename Type, bool TF>
    struct EigenBaseSelector<DynamicVector<Type, TF>>
    {
        using type = Eigen::Matrix<
            Type, 
            TF == columnVector ? Eigen::Dynamic : 1, 
            TF == rowVector ? Eigen::Dynamic : 1
        >;
    };

    template <typename Type, bool TF = defaultTransposeFlag>
    struct DynamicVector
    :   EigenBase<DynamicVector<Type, TF>>
    ,   Vector<DynamicVector<Type, TF>, TF>
    {
        typedef EigenBase<DynamicVector<Type, TF>> Base;
        typedef Type ElementType;

        DynamicVector()
        {            
        }

        explicit DynamicVector(size_t M)
        :   Base(M)
        {        
        }

        explicit DynamicVector(size_t M, Type const& val)
        :   Base(M)
        {        
            this->setConstant(val);
        }

        explicit DynamicVector(initializer_list<Type> list)
        :   Base(list.size())
        {
            assign(*this, list);
        }

        template <typename T>
        DynamicVector(Eigen::MatrixBase<T> const& rhs)
        :   Base(rhs)
        {
        }

        template <typename T>
        DynamicVector& operator=(Eigen::MatrixBase<T> const& rhs)
        {
            Base::operator=(rhs);
            return *this;
        }

        DynamicVector& operator=(ElementType const& rhs)
        {
            this->setConstant(rhs);
            return *this;
        }

        template <typename T1>
        DynamicVector& operator=(initializer_list<T1> list)
        {
            assign(*this, list);
            return *this;
        }
    };

    template <typename VT, bool TF>
    struct Rand<DynamicVector<VT, TF>>
    {
        decltype(auto) generate(size_t N)
        {
            return DynamicVector<VT, TF>::Random(N);
        }
    };

    // GET RID OF
    /*
    template <typename Type, bool TF>
    struct EigenType<DynamicVector<Type, TF>>
    {
        typedef typename DynamicVector<Type, TF>::Base type;
    };
    */
}