#pragma once

#include <tmpc/SizeT.hpp>

#include <blaze/Math.h>


namespace tmpc
{
    template <typename Real>
    class CentralDifferenceDerivative
    {
    public:
        /// @brief Constructor
        ///
        /// @param nx dimensionality of the function's domain
        /// @param ny dimensionality of the function's codomain
        CentralDifferenceDerivative(size_t nx, size_t ny)
        :   nX_(nx)
        ,   nY_(ny)
        ,   x_(nx)
        ,   yPlus_(ny)
        ,   yMinus_(ny)
        ,   J_(ny, nx)
        {
        }


        /// @Calculate numeric Jacobian df(x)/dx at point x.
        ///
        /// @param f functor object representing a map from R^NX to R^NY
        /// @param x point at which to calculate the Jacobian
        /// @param delta vector defining finite difference step in each direction
        template <typename Func, typename VT1, typename VT2>
        auto const& operator()(Func const& f, blaze::Vector<VT1, blaze::columnVector> const& x, blaze::Vector<VT2, blaze::columnVector> const& delta)
        {
            x_ = x;

            for (size_t i = 0; i < nX_; ++i)
            {
                x_[i] = (~x)[i] - (~delta)[i];
                yMinus_ = f(x_);

                x_[i] = (~x)[i] + (~delta)[i];
                yPlus_ = f(x_);

                x_[i] = (~x)[i];

                column(J_, i) = (yPlus_ - yMinus_) / (2 * (~delta)[i]);
            }

            return J_;
        }


    private:
        size_t const nX_;
        size_t const nY_;

        blaze::DynamicVector<Real, blaze::columnVector> mutable x_;
        blaze::DynamicVector<Real, blaze::columnVector> mutable yPlus_;
        blaze::DynamicVector<Real, blaze::columnVector> mutable yMinus_;
        blaze::DynamicMatrix<Real> mutable J_;
    };
}