#pragma once

#include "StorageOrder.hpp"
#include "EigenType.hpp"
#include "InitializerList.hpp"
#include "MatrixAssign.hpp"
#include "EigenBase.hpp"

#include <Eigen/Dense>

#include <type_traits>

namespace tmpc 
{

    template <typename MT, bool SO>
    struct Matrix : EigenBase<MT>
    {
        MT& operator~()
        {
            return static_cast<MT&>(*this);
        }

        MT const& operator~() const
        {
            return static_cast<MT const&>(*this);
        }

        typedef EigenBase<MT> Base;
        typedef typename Base::Scalar ElementType;

    protected:
        Matrix(Base&& rhs)
        :   Base(rhs)
        {            
        }
        
        Matrix(size_t M, size_t N)
        :   Base(M, N)
        {        
        }

        Matrix(size_t M, size_t N, ElementType const& val)
        :   Base(M, N)
        {
            this->setConstant(val);
        }

        Matrix(initializer_list<initializer_list<ElementType>> list)
        :   Base(list.size(), determineColumns(list))
        {
            assign(*this, list);
        }

        template <typename T>
        Matrix(Eigen::MatrixBase<T> const& rhs)
        :   Base(rhs)
        {        
        }
    };

    template <typename MT>
    size_t rows(Eigen::MatrixBase<MT> const& m)
    {
        return m.rows();
    }

    template <typename MT>
    size_t columns(Eigen::MatrixBase<MT> const& m)
    {
        return m.cols();
    }
}