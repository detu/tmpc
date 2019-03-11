#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/Math.hpp>

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
        ,   covariance_(n)
        ,   sample_(n)
        {
            covariance(blaze::IdentityMatrix<Real>(n));
        }


        /// @brief Creates a multinormal distribution with given mean and covariance.
        ///
        /// @param n number of dimensions
        template <typename VT, typename MT, bool SO>
        MultivariateNormalDistribution(blaze::Vector<VT, blaze::columnVector> const& mean, blaze::Matrix<MT, SO> const& cov)
        {
            auto const n = size(mean);

            if (rows(cov) != n || columns(cov) != n)
                throw std::invalid_argument("Inconsistent dimensions of mean and covariance in MultivariateNormalDistribution ctor");

            mean_ = mean;
            sample_.resize(n);
            covariance_.resize(n);
            covariance(cov);
        }


        template <typename Generator>
        auto const& operator()(Generator& g)
        {
            // Generate a vector of I.I.D. normal values
            for (auto& x : sample_)
                x = uninormal_(g);

            // Multiply by Cholesky factor of the covariance
            sample_ = covarianceFactor_ * sample_;

            // Add mean
            sample_ += mean_;

            return sample_;
        }


        /// @brief Number of the distribution dimensions.
        auto const dimension() const
        {
            return size(sample_);
        }


        /// @brief Set covariance matrix
        template <typename MT>
        void covariance(blaze::SymmetricMatrix<MT> const& val)
        {
            noresize(covariance_) = val;
            llh(covariance_, covarianceFactor_);
        }


        /// @brief Set covariance matrix
        template <typename MT, bool SO>
        void covariance(blaze::Matrix<MT, SO> const& val)
        {
            noresize(covariance_) = val;
            llh(covariance_, covarianceFactor_);
        }


        /// @brief Get covariance matrix
        auto const& covariance() const
        {
            return covariance_;
        }


        /// @brief Set mean
        template <typename VT>
        void mean(blaze::Vector<VT, blaze::columnVector> const& val)
        {
            noresize(mean_) = val;
        }


        /// @brief Get mean
        auto const& mean() const
        {
            return mean_;
        }


    private:
        blaze::DynamicVector<Real, blaze::columnVector> mean_;
        blaze::SymmetricMatrix<blaze::DynamicMatrix<Real>> covariance_;

        // Cholesky factor of the covariance matrix
        blaze::LowerMatrix<blaze::DynamicMatrix<Real>> covarianceFactor_;
        blaze::DynamicVector<Real> sample_;

        std::normal_distribution<Real> uninormal_;
    };
}