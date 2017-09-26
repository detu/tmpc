#pragma once

#include <tmpc/matrix/TransposeFlag.hpp>
#include "IsVector.hpp"
#include "IsEigenMatrix.hpp"

#include <type_traits>

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor 
{

    template <typename VT, TransposeFlag TF>
    struct Vector
    {
        VT& operator~()
        {
            return static_cast<VT&>(*this);
        }

        VT const& operator~() const
        {
            return static_cast<VT const&>(*this);
        }
    };

    

    template <typename VT, TransposeFlag TF>
    inline size_t size(Vector<VT, TF> const& v)
    {
        return (~v).size();
    }

}

namespace Eigen
{
    template <typename VT>
    inline std::enable_if_t<::tmpc::eigen_adaptor::IsEigenMatrix<VT>::value, size_t> size(VT const& v)
    {
        return v.size();
    }
}