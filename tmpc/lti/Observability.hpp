#pragma once

#include <tmpc/math/Rank.hpp>

#include <blaze/Math.h>

#include <boost/throw_exception.hpp>

#include <stdexcept>


namespace tmpc
{
    /// @brief Computes observability matrix a continuous- or discrete-time LTI system.
    ///
    /// @param A \a A matrix of the system
    /// @param C \a C matrix of the system
    /// @param O the resulting system observability matrix.
    template <typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3>
    inline void observabilityMatrix(blaze::Matrix<MT1, SO1> const& A, blaze::Matrix<MT2, SO2> const& C, blaze::Matrix<MT3, SO3>& O)
    {
        using ET = blaze::ElementType_t<MT2>;
        auto const ny = rows(C);
        auto const nx = columns(C);

        if (rows(A) != nx || columns(A) != nx)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Inconsistent matrix size"));

        resize(~O, ny * nx, nx);

        for (size_t i = 0; i < nx; ++i)
        {
            if (i == 0)
                submatrix(O, ny * i, 0, ny, nx) = ~C;
            else
                submatrix(O, ny * i, 0, ny, nx) = submatrix(O, ny * (i - 1), 0, ny, nx) * ~A;
        }
    }


    /// @brief Check if a continuous- or discrete-time LTI system is observable.
    ///
    /// @param A \a A matrix of the system
    /// @param C \a C matrix of the system
    /// @return true iff the system is observable.
    template <typename MT1, bool SO1, typename MT2, bool SO2>
    inline bool isObservable(blaze::Matrix<MT1, SO1> const& A, blaze::Matrix<MT2, SO2> const& C)
    {
        using ET = blaze::ElementType_t<MT2>;

        auto const nx = columns(C);

        if (rows(A) != nx || columns(A) != nx)
            BOOST_THROW_EXCEPTION(std::invalid_argument("Inconsistent matrix size"));

        blaze::DynamicMatrix<ET> O;
        observabilityMatrix(A, C, O);

        return rank(O) == nx;
    }
}