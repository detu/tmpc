#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/Math.hpp>

#include <stdexcept>


namespace tmpc
{
    /// @brief Defines a group of vector elements.
    ///
    /// Example:
    /// 
    /// auto constexpr stateQ = ElementGroup<0, 1> {};
    /// auto constexpr stateTheta = ElementGroup<stateQ.end, 1> {};
    /// auto constexpr stateV = ElementGroup<stateTheta.end, 1> {};
    /// auto constexpr stateOmega = ElementGroup<stateV.end, 1> {};
    ///
    template <size_t From, size_t Len>
    struct ElementGroup final
    {
        /// @brief Index of first element
        static size_t constexpr begin = From;

        /// @brief Number of elements
        static size_t constexpr size = Len;

        /// @brief Index of next-to-last element
        static size_t constexpr end = From + Len;
    };
    
    
    /// @brief Get submatrix corresponding to a pair of element groups.
    ///
    /// @param m the matrix from which the submatrix is taken
    /// @param c1 element group defining rows of the returned submatrix
    /// @param c2 element group defining columns of the returned submatrix
    /// @return submatrix of \a m with rows defined by \a c1 and columns defined by \a c2
    template <typename MT, bool SO, size_t From1, size_t Len1, size_t From2, size_t Len2>
    inline decltype(auto) submatrix(blaze::Matrix<MT, SO>& m, ElementGroup<From1, Len1> c1, ElementGroup<From2, Len2> c2)
    {
        return blaze::submatrix<From1, From2, Len1, Len2>(m);
    }


    /// @brief Get const submatrix corresponding to a pair of vector element groups.
    ///
    /// @param m the matrix from which the submatrix is taken
    /// @param c1 element group defining rows of the returned submatrix
    /// @param c2 element group defining columns of the returned submatrix
    /// @return submatrix of \a m with rows defined by \a c1 and columns defined by \a c2
    template <typename MT, bool SO, size_t From1, size_t Len1, size_t From2, size_t Len2>
    inline decltype(auto) submatrix(blaze::Matrix<MT, SO> const& m, ElementGroup<From1, Len1>, ElementGroup<From2, Len2>)
    {
        return blaze::submatrix<From1, From2, Len1, Len2>(m);
    }


    /// @brief Get subvector corresponding to an element group.
    ///
    /// @param v the vector from which the subvector is taken
    /// @param c element group defining elements of the returned subvector
    /// @return subvector of \a v with elements defined by \a c
    template <typename VT, bool TF, size_t From, size_t Len>
    inline decltype(auto) subvector(blaze::Vector<VT, TF>& v, ElementGroup<From, Len> c)
    {
        return blaze::subvector<From, Len>(v);
    }


    /// @brief Get const subvector corresponding to an element group.
    ///
    /// @param v the vector from which the subvector is taken
    /// @param c element group defining elements of the returned subvector
    /// @return subvector of \a v with elements defined by \a c
    template <typename VT, bool TF, size_t From, size_t Len>
    inline decltype(auto) subvector(blaze::Vector<VT, TF> const& v, ElementGroup<From, Len>)
    {
        return blaze::subvector<From, Len>(v);
    }
}