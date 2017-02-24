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

        typedef EigenBase<MT> OurEigenBase;
        typedef typename OurEigenBase::Scalar ElementType;

    protected:
        Matrix()
        {            
        }

        Matrix(OurEigenBase&& rhs)
        :   OurEigenBase(rhs)
        {            
        }

        Matrix(OurEigenBase const& rhs)
        :   OurEigenBase(rhs)
        {            
        }
        
        Matrix(size_t M, size_t N)
        :   OurEigenBase(M, N)
        {        
        }

        Matrix(size_t M, size_t N, ElementType const& val)
        :   OurEigenBase(M, N)
        {
            this->setConstant(val);
        }

        Matrix(initializer_list<initializer_list<ElementType>> list)
        :   OurEigenBase(list.size(), determineColumns(list))
        {
            assign(*this, list);
        }

        template <typename T>
        Matrix(Eigen::MatrixBase<T> const& rhs)
        :   OurEigenBase(rhs)
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