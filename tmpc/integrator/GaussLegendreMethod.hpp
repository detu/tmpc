#pragma once

#include <tmpc/Exception.hpp>


namespace tmpc
{
    /// @brief Gauss-Legendre method with a specified number of steps.
    struct GaussLegendreMethod
    {
        explicit GaussLegendreMethod(size_t n_stages)
        :   stages_(n_stages)
        {
        }


        template <typename MT, bool SO, typename VT1, typename VT2>
        void butcherTableau(
            blaze::Matrix<MT, SO>& A,
            blaze::Vector<VT1, blaze::rowVector>& b,
            blaze::Vector<VT2, blaze::columnVector>& c) const
        {
            switch (stages_)
            {
                case 2:
                    // The tables for the 2- and 3-point Gauss-Legendre method are taken from 
                    // https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Gauss%E2%80%93Legendre_methods
                    ~A = {
                        {1./4., 1./4. - sqrt(3.)/6.},
                        {1./4. + sqrt(3.)/6., 1./4.}
                    };
                    ~b = {1./2., 1./2.};
                    ~c = {1./2. - sqrt(3.)/6., 1./2. + sqrt(3.)/6.};
                    
                    break;

                case 3:
                    ~A = {
                        {5./36., 2./9. - sqrt(15.)/15., 5./36. - sqrt(15.)/30.},
                        {5./36. + sqrt(15.)/24., 2./9., 5./36. - sqrt(15.)/24.},
                        {5./36. + sqrt(15.)/30., 2./9. + sqrt(15.)/15., 5./36.},
                    };
                    ~b = {5./18., 4./9., 5./18.};
                    ~c = {1./2. - sqrt(15.)/10., 1./2., 1./2. + sqrt(15.)/10.};

                    break;

                default:
                    TMPC_THROW_EXCEPTION(std::invalid_argument("Unsupported number of Gauss-Lenendre method stages"));
            }
        }


        size_t stages() const
        {
            return stages_;
        }


    private:
        size_t const stages_;
    };
}