#pragma once

#include "InitializerList.hpp"

#include "Eigen.hpp"

namespace tmpc :: eigen_adaptor
{
    template <typename T, typename T1>
    void assign(Eigen::DenseBase<T>& matrix, initializer_list<initializer_list<T1>> list)
    {
        matrix.resize(list.size(), determineColumns(list));

        size_t i = 0;
        for (auto const& row : list) 
        {
            size_t j = 0;

            for (auto const& val : row)
                matrix(i, j++) = val;

            for (; j < matrix.cols(); ++j)
                matrix(i, j) = T1();

            ++i;
        }
    }
}