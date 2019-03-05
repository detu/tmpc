#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/Math.hpp>

#include <stdexcept>


namespace tmpc
{
    /// @brief Defines ordering of vector components.
    ///
    /// Example:
    /// 
    /// using StateQ = VectorComponent<0, 1>;
    /// using StateTheta = VectorComponent<StateQ::end, 1>;
    /// using StateV = VectorComponent<StateTheta::end, 1>;
    /// using StateOmega = VectorComponent<StateV::end, 1>;
    ///
    template <size_t From, size_t Len>
    struct VectorComponent
    {
        static size_t constexpr begin = From;
        static size_t constexpr size = Len;
        static size_t constexpr end = From + Len;
    };


    namespace detail
    {
        template <typename C>
        struct VectorComponentImpl;


        template <size_t From, size_t Len>
        struct VectorComponentImpl<VectorComponent<From, Len>>
        {
            template <typename VT, bool TF>
            static decltype(auto) component(blaze::Vector<VT, TF>& v)
            {
                return blaze::subvector<From, Len>(v);
            }


            template <typename VT, bool TF>
            static decltype(auto) component(blaze::Vector<VT, TF> const& v)
            {
                return blaze::subvector<From, Len>(v);
            }
        };


        template <typename C1, typename C2>
        struct MatrixComponentImpl;


        template <size_t From1, size_t Len1, size_t From2, size_t Len2>
        struct MatrixComponentImpl<VectorComponent<From1, Len1>, VectorComponent<From2, Len2>>
        {
            template <typename MT, bool SO>
            static decltype(auto) component(blaze::Matrix<MT, SO>& m)
            {
                return blaze::submatrix<From1, From2, Len1, Len2>(m);
            }


            template <typename MT, bool SO>
            static decltype(auto) component(blaze::Matrix<MT, SO> const& m)
            {
                return blaze::submatrix<From1, From2, Len1, Len2>(m);
            }
        };
    }


    template <typename C1, typename C2, typename MT, bool SO>
    inline decltype(auto) component(blaze::Matrix<MT, SO>& m)
    {
        return detail::MatrixComponentImpl<C1, C2>::component(m);
    }


    template <typename C1, typename C2, typename MT, bool SO>
    inline decltype(auto) component(blaze::Matrix<MT, SO> const& m)
    {
        return detail::MatrixComponentImpl<C1, C2>::component(m);
    }


    template <typename C, typename VT, bool TF>
    inline decltype(auto) component(blaze::Vector<VT, TF>& v)
    {
        return detail::VectorComponentImpl<C>::component(v);
    }


    template <typename C, typename VT, bool TF>
    inline decltype(auto) component(blaze::Vector<VT, TF> const& v)
    {
        return detail::VectorComponentImpl<C>::component(v);
    }
}