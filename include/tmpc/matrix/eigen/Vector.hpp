#pragma once

#include "TransposeFlag.hpp"
#include "IsVector.hpp"
#include "EigenType.hpp"

#include <type_traits>

#include <Eigen/Dense>

namespace tmpc 
{

    template <typename VT, bool TF>
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

    template <typename T>
    inline std::enable_if_t<IsVector<T>::value, size_t> size(Eigen::MatrixBase<T> const& v)
    {
        return v.size();
    }

    template <typename VT, bool TF>
    inline size_t size(Vector<VT, TF> const& v)
    {
        return (~v).size();
    }

}