#pragma once

#include "AlignmentFlag.hpp"
#include "EigenType.hpp"
#include "Matrix.hpp"

#include <Eigen/Dense>

namespace tmpc {

/*------------------------------------------------------
 *
 * Submatrices
 *
 -------------------------------------------------------*/
template <typename MT, bool AF>
struct Submatrix;

template <typename MT, bool AF>
using SubmatrixBase = Matrix<Submatrix<MT, AF>, MT::IsRowMajor ? rowMajor : columnMajor>;

template <typename MT, bool AF = unaligned>
struct Submatrix
:   SubmatrixBase<MT, AF>
{
    typedef SubmatrixBase<MT, AF> Base;

    Submatrix(EigenBase<MT>&& rhs)
    :   Base(rhs)
    {        
    }

    Submatrix(EigenBase<MT> const& rhs)
    :   Base(rhs)
    {        
    }

    Submatrix(Submatrix const& rhs) = default;
};

template <typename MT>
decltype(auto) submatrix(Eigen::MatrixBase<MT>& matrix, size_t row, size_t column, size_t m, size_t n)
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