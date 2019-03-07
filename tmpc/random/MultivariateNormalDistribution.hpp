#pragma once

#include <tmpc/SizeT.hpp>

#include <blaze/Math.h>

#include <algorithm>
#include <random>
#include <stdexcept>


namespace tmpc
{
    /// @brief Generates random numbers according to the Multivariate Normal (or Gaussian) random number distribution.
    ///
    /// https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    ///
    template <typename Real>
    class MultivariateNormalDistribution
    {
    public:
        /// @brief Creates a multinormal distribution with zero mean and identity covariance.
        ///
        /// @param n number of dimensions
        MultivariateNormalDistribution(size_t n)
        :   mean_(blaze::ZeroVector<Real>(n))
        ,   covariance_(blaze::IdentityMatrix<Real>(n))
        ,   sample_(n)
        {
        }


        /// @brief Creates a multinormal distribution with given mean and covariance.
        ///
        /// @param n number of dimensions
        template <typename VT, typename MT, bool SO>
        MultivariateNormalDistribution(blaze::Vector<VT, blaze::columnVector> const& mean, blaze::Matrix<MT, SO> const& covariance)
        {
            auto const n = size(mean);

            if (rows(covariance) != n || columns(covariance) != n)
                throw std::invalid_argument("Inconsistent dimensions of mean and covariance in MultivariateNormalDistribution ctor");

            mean_ = mean;
            covariance_ = covariance;
            sample_.resize(n);
        }


        template <typename Generator>
        auto const& operator()(Generator& g)
        {
            // Generate a vector of I.I.D. normal values
            for (auto& x : sample_)
                x = uninormal_(g);

            // ...
            return sample_;
        }


        /// @brief Number of the distribution dimensions.
        auto const dimension() const
        {
            return size(sample_);
        }


    private:
        blaze::DynamicVector<Real, blaze::columnVector> mean_;
        blaze::DynamicMatrix<Real> covariance_;
        blaze::DynamicVector<Real> sample_;

        std::normal_distribution<Real> uninormal_;
    };
}