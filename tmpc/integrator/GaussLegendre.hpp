#pragma once

#include <tmpc/integrator/ButcherTableau.hpp>

#include <stdexcept>


namespace tmpc
{
    /// @brief Butcher tableau for Gauss-Legendre quadrature with a given number of steps.
    template <typename Real>
    inline ButcherTableau<Real> gaussLegendre(size_t n_steps)
    {
        using Tableau = ButcherTableau<Real>;

        if (n_steps != 2)
            throw std::invalid_argument(std::string(__func__) + ": only 2-point Gauss-Lenendre quadrature is currently supported");

        // The table for 2-point method is taken from 
        // https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Gauss%E2%80%93Legendre_methods
        return Tableau(
            typename Tableau::A_type {
                {1./4., 1./4. - sqrt(3.)/6.},
                {1./4. + sqrt(3.)/6., 1./4.}
            },
            typename Tableau::b_type {1./2., 1./2.},
            typename Tableau::c_type {1./2. - sqrt(3.)/6., 1./2. + sqrt(3.)/6.}
        );
    }
}