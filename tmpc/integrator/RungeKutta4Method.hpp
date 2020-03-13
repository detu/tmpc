#pragma once

#include <blaze/Math.h>


namespace tmpc
{
    /// @brief Runge-Kutta method of order 4.
    struct RungeKutta4Method
    {
        template <typename MT, bool SO, typename VT1, typename VT2>
        void butcherTableau(
            blaze::Matrix<MT, SO>& A,
            blaze::Vector<VT1, blaze::rowVector>& b,
            blaze::Vector<VT2, blaze::columnVector>& c) const
        {
            ~A = {
				{0., 0., 0., 0.},
				{1./2., 0., 0., 0.},
				{0., 1./2., 0., 0.},
				{0., 0., 1., 0.}
			};

            ~b = {1./6., 1./3., 1./3., 1./6.};
            ~c = {0., 1./2., 1./2., 1.};
        }


        size_t stages() const
        {
            return 4;
        }
    };
}