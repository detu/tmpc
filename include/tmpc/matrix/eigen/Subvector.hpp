#pragma once

#include "AlignmentFlag.hpp"
#include "IsVector.hpp"
#include "EigenType.hpp"

#include "Eigen.hpp"

#include <type_traits>

namespace tmpc {

/*------------------------------------------------------
 *
 * Subvectors
 *
 -------------------------------------------------------*/
template <typename VT>
using UnalignedSubvectorBase = Eigen::Block<
        typename EigenType<VT>::type, 
        VT::RowsAtCompileTime == 1 ? 1 : Eigen::Dynamic, 
        VT::ColsAtCompileTime == 1 ? 1 : Eigen::Dynamic
    >;

template <typename VT, bool AF = unaligned, typename Enable = void>
struct Subvector;

template <typename VT>
struct Subvector<VT, unaligned, std::enable_if_t<IsVector<VT>::value>> 
:   UnalignedSubvectorBase<VT>
{
    typedef UnalignedSubvectorBase<VT> Base;
    typedef typename Base::Scalar ElementType;

    Subvector(Base&& rhs)
    :   Base(rhs)
    {        
    }

    Subvector& operator=(ElementType const& rhs)
    {
        this->setConstant(rhs);
        return *this;
    }

    template <typename T>
    Subvector& operator=(Eigen::MatrixBase<T> const& rhs)
    {
        Base::operator=(rhs);
        return *this;
    }
};

template <typename VT>
Subvector<VT, unaligned> subvector(VT& vector, size_t index, size_t size)
{
    return vector.segment(index, size);
}

}