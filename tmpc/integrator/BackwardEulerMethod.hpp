#pragma once

#include <blaze/Math.h>


namespace tmpc
{
    /// @brief Backward Euler method.
    struct BackwardEulerMethod
    {
        template <typename MT, bool SO, typename VT1, typename VT2>
        void butcherTableau(
            blaze::Matrix<MT, SO>& A,
            blaze::Vector<VT1, blaze::rowVector>& b,
            blaze::Vector<VT2, blaze::columnVector>& c) const
        {
            ~A = {{1.}};
            ~b = {1.};
            ~c = {1.};
        }


        size_t stages() const
        {
            return 1;
        }
    };
}