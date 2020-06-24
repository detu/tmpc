#pragma once

#include <tmpc/math/Rank.hpp>

#include <blaze/Math.h>

#include <boost/throw_exception.hpp>

#include <stdexcept>


namespace tmpc
{
    /// @brief Computes controllability matrix a continuous- or discrete-time LTI system.
    ///
    /// @param A \a A matrix of the system
    /// @param B \a B matrix of the system
    /// @param C the resulting system controllability matrix.
    template <typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3>
    inline auto controllabilityMatrix(blaze::Matrix<MT1, SO1> const& A, blaze::Matrix<MT2, SO2> const& B, blaze::Matrix<MT3, SO3>& C)
    {
        using ET = blaze::ElementType_t<MT1>;
        auto const nx = rows(B);
        auto const nu = columns(B);

        if (rows(A) != nx || columns(A) != nx)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Inconsistent matrix size"));

        resize(~C, nx, nx * nu);

        for (size_t i = 0; i < nx; ++i)
        {
            if (i == 0)
                submatrix(C, 0, nu * i, nx, nu) = ~B;
            else
                submatrix(C, 0, nu * i, nx, nu) = ~A * submatrix(C, 0, nu * (i - 1), nx, nu);
        }
    }


    /// @brief Check if a continuous- or discrete-time LTI system is controllable.
    ///
    /// @param A \a A matrix of the system
    /// @param B \a B matrix of the system
    /// @return true iff the system is controllable.
    template <typename MT1, bool SO1, typename MT2, bool SO2>
    inline bool isControllable(blaze::Matrix<MT1, SO1> const& A, blaze::Matrix<MT2, SO2> const& B)
    {
        using ET = blaze::ElementType_t<MT1>;

        auto const nx = rows(B);

        if (rows(A) != nx || columns(A) != nx)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Inconsistent matrix size"));

        blaze::DynamicMatrix<ET> C;
        controllabilityMatrix(A, B, C);

        return rank(C) == nx;
    }
}