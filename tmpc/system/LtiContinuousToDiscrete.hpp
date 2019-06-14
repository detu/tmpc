#pragma once

#include <tmpc/Math.hpp>


namespace tmpc
{
    /// @brief Convert continuous time linear state space dynamics to discrete time.
    ///
    /// For the formulas, look, for example, here: http://www.engr.iupui.edu/~skoskie/ECE595_f05/handouts/discretization.pdf
    /// The idea of the single matrix exponential algorithm is taken from 
    /// https://github.com/scipy/scipy/blob/f07bfba/scipy/signal/lti_conversion.py#L449
    ///
    template <typename Real>
    class LtiContinuousToDiscrete
    {
    public:
        LtiContinuousToDiscrete(size_t nx, size_t nu)
        :   nx_(nx)
        ,   nu_(nu)
        ,   work_(nx + nu, nx + nu)
        {
        }


        template <typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3, typename MT4, bool SO4>
        void operator()(blaze::Matrix<MT1, SO1> const& Ac, blaze::Matrix<MT2, SO2> const& Bc, Real time_step,
            blaze::Matrix<MT3, SO3>& Ad, blaze::Matrix<MT4, SO4>& Bd)
        {
            submatrix(work_, 0, 0, nx_, nx_) = Ac;
            submatrix(work_, 0, nx_, nx_, nu_) = Bc;
            submatrix(work_, nx_, 0, nu_, nx_ + nu_) = 0.;

            work_ = expm(time_step * work_);

            ~Ad = submatrix(work_, 0, 0, nx_, nx_);
            ~Bd = submatrix(work_, 0, nx_, nx_, nu_);
        }


    private:
        size_t nx_;
        size_t nu_;
        blaze::DynamicMatrix<Real> work_;
    };
}