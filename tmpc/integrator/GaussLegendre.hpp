#pragma once

#include <tmpc/integrator/ButcherTableau.hpp>

#include <boost/throw_exception.hpp>

#include <stdexcept>


namespace tmpc
{
    /// @brief Butcher tableau for Gauss-Legendre quadrature with a given number of steps.
    template <typename Real>
    inline ButcherTableau<Real> gaussLegendre(size_t n_steps)
    {
        using Tableau = ButcherTableau<Real>;

        switch (n_steps)
        {
            case 2:
                // The tables for the 2- and 3-point Gauss-Legendre method are taken from 
                // https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Gauss%E2%80%93Legendre_methods
                return Tableau(
                    typename Tableau::A_type {
                        {1./4., 1./4. - sqrt(3.)/6.},
                        {1./4. + sqrt(3.)/6., 1./4.}
                    },
                    typename Tableau::b_type {1./2., 1./2.},
                    typename Tableau::c_type {1./2. - sqrt(3.)/6., 1./2. + sqrt(3.)/6.}
                );
            case 3:
                return Tableau(
                    typename Tableau::A_type {
                        {5./36., 2./9. - sqrt(15.)/15., 5./36. - sqrt(15.)/30.},
                        {5./36. + sqrt(15.)/24., 2./9., 5./36. - sqrt(15.)/24.},
                        {5./36. + sqrt(15.)/30., 2./9. + sqrt(15.)/15., 5./36.},
                    },
                    typename Tableau::b_type {5./18., 4./9., 5./18.},
                    typename Tableau::c_type {1./2. - sqrt(15.)/10., 1./2., 1./2. + sqrt(15.)/10.}
                );
            default:
                BOOST_THROW_EXCEPTION(std::invalid_argument("Only 2-point Gauss-Lenendre quadrature is currently supported"));
        }
    }
}