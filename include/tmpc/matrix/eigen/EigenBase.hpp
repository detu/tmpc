#pragma once

#include "Forward.hpp"

#include "Eigen.hpp"

namespace tmpc
{
    template <typename T>
    struct EigenBaseSelector;

    template <typename T>
    using EigenBase = typename EigenBaseSelector<T>::type;

    template <typename Type, size_t M, size_t N, bool SO>
    struct EigenBaseSelector<StaticMatrix<Type, M, N, SO>>
    {
        static bool constexpr storageOrder = 
            (M > 1 && N == 1) 
            ? Eigen::ColMajor 
            : (SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor);

        using type = Eigen::Matrix<Type, M, N, storageOrder>;
    };

    template <typename Type, bool SO>
    struct EigenBaseSelector<DynamicMatrix<Type, SO>>
    {
        using type = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic, 
            SO == rowMajor ? Eigen::RowMajor : Eigen::ColMajor>;
    };

    template <typename MT, bool AF>
    struct EigenBaseSelector<Submatrix<MT, AF>>
    {
        using type = Eigen::Block<EigenBase<MT>, Eigen::Dynamic, Eigen::Dynamic>;
    };
}