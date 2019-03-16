#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/Math.hpp>


namespace tmpc
{
    /// @brief Kalman filter implementation.
    ///
    template <typename Real>
    class KalmanFilter
    {
    public:
        using Matrix = blaze::DynamicMatrix<Real>;
        using SymmetricMatrix = blaze::SymmetricMatrix<Matrix>;


        KalmanFilter(size_t nx, size_t ny)
        :   nx_(nx)
        ,   ny_(ny)
        ,   stateEstimate_(nx, Real {})
        ,   stateCovariance_(nx)
        ,   processNoiseCovariance_(nx)
        ,   measurementNoiseCovariance_(ny)
        ,   S_(ny)
        ,   S_chol_(ny)
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
        SymmetricMatrix const& stateCovariance() const
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
        SymmetricMatrix const& processNoiseCovariance() const
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
        SymmetricMatrix const& measurementNoiseCovariance() const
        {
            return measurementNoiseCovariance_;
        }


        /// @brief Set measurement noise covariance
        template <typename MT, bool SO>
        void measurementNoiseCovariance(blaze::Matrix<MT, SO> const& val)
        {
            noresize(measurementNoiseCovariance_) = val;
        }


        /// @brief Update state estimate based on next state value and sensitivities.
        ///
        /// @param x_next next state value
        /// @param A sensitivity matrix, A = d(x_next)/dx
        template <typename VT, typename MT, bool SO>
        void predict(blaze::Vector<VT, blaze::columnVector> const& x_next, blaze::Matrix<MT, SO> const& A)
        {
            stateEstimate_ = x_next;
            stateCovariance_ = ~A * stateCovariance_ * trans(~A) + processNoiseCovariance_;
        }


        /// @brief Update state estimate based on a linear model and control input.
        ///
        /// @param A linear model matrix A
        /// @param B linear model matrix B
        /// @param u control input
        template <typename MT1, bool SO1, typename MT2, bool SO2, typename VT>
        void predict(blaze::Matrix<MT1, SO1> const& A, blaze::Matrix<MT2, SO2> const& B, blaze::Vector<VT, blaze::columnVector> const& u)
        {
            stateEstimate_ = ~A * stateEstimate_ + ~B * ~u;
            stateCovariance_ = ~A * stateCovariance_ * trans(~A) + processNoiseCovariance_;
        }


        /// @brief Update state estimate based on measurement residual and sensitivities.
        ///
        /// @param y measurement residual, the difference between measured and predicted system output.
        /// @param C output sensitivity matrix, C = d(y_pred)/dx, where y_pred is the predicted system output for state x.
        template <typename VT, typename MT, bool SO>
        void update(blaze::Vector<VT, blaze::columnVector> const& y, blaze::Matrix<MT, SO> const& C)
        {
            // TODO: performance of this code can definitely be improved by avoiding temporary subexpressions
            // and using proper LAPACK routines.
            S_ = measurementNoiseCovariance_ + ~C * stateCovariance_ * trans(~C);
            llh(S_, S_chol_);

            // Use the Cholesky decomposition of S to multiply with inv(S).
            // K_ = stateCovariance_ * trans(~C) * inv(S_);
            K_ = stateCovariance_ * trans(~C);
            K_ *= inv(trans(S_chol_));
            K_ *= inv(S_chol_);

            // This results in the following message:
            //  ** On entry to DTRSM parameter number  9 had an illegal value
            // trsm(S_chol_, K_, CblasLeft, CblasLower, Real {1.});
            // trsm(S_chol_, K_, CblasRight, CblasLower, Real {1.});

            stateEstimate_ += K_ * ~y;

            // Use the Cholesky decomposition of S to enforce symmetry by putting the expression in the form A*A^T.
            stateCovariance_ -= K_ * S_chol_ * trans(K_ * S_chol_);
        }


        /// @brief Number of states
        auto const nx() const
        {
            return nx_;
        }


        /// @brief Number of outputs
        auto const ny() const
        {
            return ny_;
        }


        /// @brief Kalman gain on the last update step
        Matrix const& gain() const
        {
            return K_;
        }


        /// @brief Measurement prediction covariance on the last update step
        SymmetricMatrix const& measurementPredictionCovariance() const
        {
            return S_;
        }

        
    private:
        size_t const nx_;
        size_t const ny_;

        blaze::DynamicVector<Real, blaze::columnVector> stateEstimate_;
        SymmetricMatrix stateCovariance_;

        SymmetricMatrix processNoiseCovariance_;
        SymmetricMatrix measurementNoiseCovariance_;

        SymmetricMatrix S_;

        // Cholesky decomposition of S_;
        blaze::LowerMatrix<Matrix> S_chol_;

        // Kalman gain
        Matrix K_;
    };
}