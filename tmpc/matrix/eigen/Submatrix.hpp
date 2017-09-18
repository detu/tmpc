#pragma once

#include "AlignmentFlag.hpp"
#include "EigenType.hpp"
#include "Matrix.hpp"
#include "EigenBase.hpp"

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor 
{

    /*------------------------------------------------------
    *
    * Submatrices
    *
    -------------------------------------------------------*/
    template <typename MT, bool AF>
    struct Submatrix;
	
	template <typename MT, bool AF>
    struct EigenBaseSelector<Submatrix<MT, AF>>
    {
        using type = Eigen::Block<EigenBase<MT>, Eigen::Dynamic, Eigen::Dynamic>;
    };

    template <typename MT, bool AF>
    using SubmatrixBase = Matrix<Submatrix<MT, AF>, MT::IsRowMajor ? rowMajor : columnMajor>;

    template <typename MT, bool AF = unaligned>
    struct Submatrix
    :   SubmatrixBase<MT, AF>
    ,   EigenBase<Submatrix<MT, AF>>
    {
        using Base = SubmatrixBase<MT, AF>;
        using OurEigenBase = EigenBase<Submatrix<MT, AF>>;
        using ElementType = typename OurEigenBase::Scalar;

        Submatrix(OurEigenBase&& rhs)
        :   OurEigenBase(std::move(rhs))
        {        
        }

        Submatrix(OurEigenBase const& rhs)
        :   OurEigenBase(rhs)
        {        
        }

        Submatrix(Submatrix const& rhs) = default;

        Submatrix& operator=(ElementType const& rhs)
        {
            this->setConstant(rhs);
            return *this;
        }

        template <typename T>
        Submatrix& operator=(Eigen::MatrixBase<T> const& rhs)
        {
            OurEigenBase::operator=(rhs);
            return *this;
        }
    };

    template <typename MT>
    Submatrix<MT, unaligned> submatrix(Eigen::MatrixBase<MT>& matrix, size_t row, size_t column, size_t m, size_t n)
    {
        return matrix.block(row, column, m, n);
    }

    template <typename MT>
    Submatrix<MT const, unaligned> submatrix(Eigen::MatrixBase<MT> const& matrix, size_t row, size_t column, size_t m, size_t n)
    {
        return matrix.block(row, column, m, n);
    }

    /*
    template <typename MT, bool AF = unaligned>
    struct Submatrix;

    template <typename MT>
    struct Submatrix<MT, unaligned> : 
        Eigen::Block<EigenType<MT>::type, Eigen::Dynamic, Eigen::Dynamic>
    {
        typedef Eigen::Block<EigenType<MT>::type, Eigen::Dynamic, Eigen::Dynamic> Base;
    };

    template <typename MT>
    Submatrix<typename EigenType<MT>::type::Derived> submatrix(MT& matrix, size_t row, size_t column, size_t m, size_t n)
    {
        return Submatrix<typename EigenType<MT>::type::Derived>(matrix, row, column, m, n);
    }

    template <typename MT>
    Submatrix<typename EigenType<MT const>::type::Derived> submatrix(MT const& matrix, size_t row, size_t column, size_t m, size_t n)
    {
        return Submatrix<typename EigenType<MT const>::type::Derived>(matrix, row, column, m, n);
    }
    */

}