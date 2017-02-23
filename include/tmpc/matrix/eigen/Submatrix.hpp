#pragma once

#include "AlignmentFlag.hpp"
#include "EigenType.hpp"

#include <Eigen/Dense>

namespace tmpc {

/*------------------------------------------------------
 *
 * Submatrices
 *
 -------------------------------------------------------*/

template <typename MT, bool AF = unaligned>
using Submatrix = Eigen::Block<typename EigenType<MT>::type, Eigen::Dynamic, Eigen::Dynamic>;

template <typename MT>
decltype(auto) submatrix(MT& matrix, size_t row, size_t column, size_t m, size_t n)
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