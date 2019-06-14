#pragma once

#include <blaze/Math.h>

#include <unsupported/Eigen/MatrixFunctions>

#include <stdexcept>

#include <boost/throw_exception.hpp>


namespace tmpc
{
    namespace detail
    {
        template <bool SO>
        struct EigenStorageOrder;


        template <>
        struct EigenStorageOrder<blaze::columnMajor>
        {
            static auto constexpr value = Eigen::ColMajor;
        };


        template <>
        struct EigenStorageOrder<blaze::rowMajor>
        {
            static auto constexpr value = Eigen::RowMajor;
        };
    };


    /// @brief Matrix exponential.
    ///
    /// A hacky and non-optimal way of implementing it using Eigen
    /// before Blaze starts supporting it: https://bitbucket.org/blaze-lib/blaze/issues/74/matrix-exponentiation
    template <typename MT, bool SO>
    inline blaze::DynamicMatrix<typename MT::ElementType, SO> expm(blaze::Matrix<MT, SO> const& m)
    {
        auto const M = rows(m);
        if (M != columns(m))
            BOOST_THROW_EXCEPTION(std::invalid_argument("Matrix must be square"));

        using Scalar = typename MT::ElementType;
        auto constexpr EigenSO = detail::EigenStorageOrder<SO>::value;
        using EigenMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, EigenSO>;
        using EigenStride = Eigen::Stride<Eigen::Dynamic, 1>;

        blaze::DynamicMatrix<Scalar, SO> tmp = m;
        Eigen::Map<EigenMatrix, Eigen::Unaligned, EigenStride> eig_mat(data(tmp), M, M, EigenStride(spacing(tmp), 1));

        // Use Eigen to calculate matrix exponent
        eig_mat = eig_mat.exp();

        return tmp;
    }
}