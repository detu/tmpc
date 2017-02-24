#pragma once

#include "InitializerList.hpp"
#include "IsVector.hpp"

#include "Eigen.hpp"

#include <type_traits>

namespace tmpc
{
    template <typename T, typename T1>
    //std::enable_if_t<IsVector<T>::value, void> 
    void
    assign(Eigen::DenseBase<T>& v, initializer_list<T1> list)
    {
        v.resize(list.size());

        size_t i = 0;
        for (auto const& val : list)
            v[i++] = val;
    }
}