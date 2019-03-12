#pragma once

#include <tmpc/SizeT.hpp>

#include <blaze/Math.h>


namespace tmpc
{
    /// @brief Kalman filter implementation.
    ///
    template <typename Real>
    class KalmanFilter
    {
    public:
        KalmanFilter(size_t nx, size_t nu, size_t ny)
        :   nx_(nx)
        ,   nu_(nu)
        ,   ny_(ny)
        ,   stateEstimate_(nx, Real {})
        ,   stateCovariance_(nx)
        ,   A_(blaze::IdentityMatrix<Real>(nx))
        ,   B_(nx, nu)
        ,   C_(ny, nx, Real {})
        ,   processNoiseCovariance_(nx)
        ,   measurementNoiseCovariance_(ny)
        ,   y_(ny)
        ,   S_(ny)
        ,   K_(nx, ny)
        {
        }


        /// @brief Get state estimate
        auto const& stateEstimate() const
        {
            return stateEstimate_;
        }


        /// @brief Set state estimate
        template <typename VT>
        void stateEstimate(blaze::Vector<VT, blaze::columnVector> const& val)
        {
            noresize(stateEstimate_) = val;
        }


        /// @brief Get state covariance
        auto const& stateCovariance() const
        {
            return stateCovariance_;
        }


        /// @brief Set state covariance
        template <typename MT, bool SO>
        void stateCovariance(blaze::Matrix<MT, SO> const& val)
        {
            noresize(stateCovariance_) = val;
        }


        /// @brief Get process noise covariance
        auto const& processNoiseCovariance() const
        {
            return processNoiseCovariance_;
        }


        /// @brief Set process noise covariance
        template <typename MT, bool SO>
        void processNoiseCovariance(blaze::Matrix<MT, SO> const& val)
        {
            noresize(processNoiseCovariance_) = val;
        }


        /// @brief Get measurement noise covariance
        auto const& measurementNoiseCovariance() const
        {
            return measurementNoiseCovariance_;
        }


        /// @brief Set measurement noise covariance
        template <typename MT, bool SO>
        void measurementNoiseCovariance(blaze::Matrix<MT, SO> const& val)
        {
            noresize(measurementNoiseCovariance_) = val;
        }


        /// @brief Get A matrix
        auto const& A() const
        {
            return A_;
        }


        /// @brief Set A matrix
        template <typename MT, bool SO>
        void A(blaze::Matrix<MT, SO> const& val)
        {
            noresize(A_) = val;
        }


        /// @brief Get B matrix
        auto const& B() const
        {
            return B_;
        }


        /// @brief Set B matrix
        template <typename MT, bool SO>
        void B(blaze::Matrix<MT, SO> const& val)
        {
            noresize(B_) = val;
        }


        /// @brief Get C matrix
        auto const& C() const
        {
            return C_;
        }


        /// @brief Set C matrix
        template <typename MT, bool SO>
        void C(blaze::Matrix<MT, SO> const& val)
        {
            noresize(C_) = val;
        }


        template <typename VT>
        void predict(blaze::Vector<VT, blaze::columnVector> const& u)
        {
            stateEstimate_ = A_ * stateEstimate_ + B_ * u;
            stateCovariance_ = A_ * stateCovariance_ * trans(A_) + processNoiseCovariance_;
        }


        template <typename VT>
        void update(blaze::Vector<VT, blaze::columnVector> const& z)
        {
            y_ = z - C_ * stateEstimate_;
            S_ = measurementNoiseCovariance_ + C_ * stateCovariance_ * trans(C_);
            K_ = stateCovariance_ * trans(C_) * inv(S_);
            stateEstimate_ += K_ * y_;
            stateCovariance_ = (blaze::IdentityMatrix<Real>(nx_) - K_ * C_) * stateCovariance_;
        }


        /// @brief Number of states
        auto const nx() const
        {
            return nx_;
        }


        /// @brief Number of inputs
        auto const nu() const
        {
            return nu_;
        }


        /// @brief Number of outputs
        auto const ny() const
        {
            return ny_;
        }

        
    private:
        size_t const nx_;
        size_t const nu_;
        size_t const ny_;

        blaze::DynamicVector<Real, blaze::columnVector> stateEstimate_;
        blaze::SymmetricMatrix<blaze::DynamicMatrix<Real>> stateCovariance_;

        blaze::DynamicMatrix<Real> A_;
        blaze::DynamicMatrix<Real> B_;
        blaze::DynamicMatrix<Real> C_;

        blaze::SymmetricMatrix<blaze::DynamicMatrix<Real>> processNoiseCovariance_;
        blaze::SymmetricMatrix<blaze::DynamicMatrix<Real>> measurementNoiseCovariance_;

        blaze::DynamicVector<Real, blaze::columnVector> y_;
        blaze::SymmetricMatrix<blaze::DynamicMatrix<Real>> S_;
        blaze::DynamicMatrix<Real> K_;
    };
}