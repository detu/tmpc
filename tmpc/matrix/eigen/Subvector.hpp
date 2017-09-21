#pragma once

#include <tmpc/matrix/AlignmentFlag.hpp>
#include "IsVector.hpp"
#include "EigenType.hpp"
#include "IsRowVector.hpp"
#include "EigenBase.hpp"

#include <type_traits>

namespace tmpc :: eigen_adaptor 
{
    template <typename VT, AlignmentFlag AF, TransposeFlag TF>
    struct Subvector;

    template <typename VT, AlignmentFlag AF, TransposeFlag TF>
    struct EigenBaseSelector<Subvector<VT, AF, TF>>
    {
        using type = Eigen::VectorBlock<
            EigenBase<VT>, 
            Eigen::Dynamic
        >;
    };

    /*------------------------------------------------------
    *
    * Subvectors
    *
    -------------------------------------------------------*/
    template <
        typename VT, 
        AlignmentFlag AF = unaligned,
        TransposeFlag TF = IsRowVector<VT>::value ? rowVector : columnVector
    >
    struct Subvector
    :   Vector<Subvector<VT, AF, TF>, TF>
    ,   EigenBase<Subvector<VT, AF, TF>>
    {
        typedef EigenBase<Subvector<VT, AF, TF>> OurEigenBase;
        typedef typename OurEigenBase::Scalar ElementType;

        Subvector(OurEigenBase&& rhs)
        :   OurEigenBase(std::move(rhs))
        {        
        }

        Subvector(Subvector<OurEigenBase>&& rhs)
        :   OurEigenBase(std::move(rhs))
        {        
        }

        Subvector(OurEigenBase const& rhs)
        :   OurEigenBase(rhs)
        {        
        }

        Subvector(Subvector<OurEigenBase> const& rhs)
        :   OurEigenBase(rhs)
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
            OurEigenBase::operator=(rhs);
            return *this;
        }
    };

    template <typename VT>
    Subvector<VT, unaligned> subvector(Eigen::MatrixBase<VT>& v, size_t index, size_t size)
    {
        return v.segment(index, size);
    }

    template <typename VT>
    Subvector<VT const, unaligned> subvector(Eigen::MatrixBase<VT> const& v, size_t index, size_t size)
    {
        return v.segment(index, size);
    }
}