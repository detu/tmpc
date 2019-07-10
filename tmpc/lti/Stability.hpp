#pragma once

#include <blaze/Math.h>


namespace tmpc
{
    /// @brief Check if a continuous-time LTI system is globally asymptotically stable.
    ///
    /// @param A \a A matrix of the system
    /// @return true if the system is stable.
    template <typename MT, bool SO>
    inline bool isGloballyAsymptoticallyStable(blaze::DenseMatrix<MT, SO> const& A)
    {
        // Compute eigenvalues
        using ET = blaze::ElementType_t<MT>;
        blaze::StaticVector<std::complex<Real>, NX> lam;
        eigen(A, lam);

        // Check that all eigenvalues are within the left complex half-plane
        return std::find_if(begin(lam), end(lam), [] (auto l) {real(l) >= 0.}) == end(lam);
    }


    /// @brief Check if a discrete-time LTI system is globally asymptotically stable.
    ///
    /// @param A \a A matrix of the system
    /// @param ts sampling time of the system. Used to indicate that the system is discrete-time.
    /// @return true if the system is stable.
    template <typename MT, bool SO, typename Real>
    inline bool isGloballyAsymptoticallyStable(blaze::DenseMatrix<MT, SO> const& A, Real ts)
    {
        // Compute eigenvalues
        using ET = blaze::ElementType_t<MT>;
        blaze::StaticVector<std::complex<Real>, NX> lam;
        eigen(A, lam);

        // Check that all eigenvalues are within the unit disk
        return std::find_if(begin(lam), end(lam), [] (auto l) {abs(l) >= 1.}) == end(lam);
    }
}