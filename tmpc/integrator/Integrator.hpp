#pragma once

#include <tmpc/SizeT.hpp>

#include <blaze/Math.h>

#include <cmath>


namespace tmpc
{
    /// @brief CRTP base class for integrators.
    ///
    template <typename Derived>
    class Integrator
    {
    public:
        Derived const& operator~() const noexcept
        {
            return static_cast<Derived const&>(*this);
        }


        Derived& operator~() noexcept
        {
            return static_cast<Derived&>(*this);
        }


    protected:
        Integrator() = default;
        Integrator(Integrator const&) = default;
        Integrator(Integrator&&) = default;
        Integrator& operator=(Integrator const&) = default;
        Integrator& operator=(Integrator&&) = default;
        ~Integrator() = default;
    };


    template <typename I, typename DE, typename Real, typename VT1, typename VT2, typename VT3>
    inline void integrate(
        Integrator<I> const& integrator, 
        DE const& de,
        Real t0, Real h, Real h_max,
        blaze::Vector<VT1, blaze::columnVector> const& x0, 
        blaze::Vector<VT2, blaze::columnVector> const& u,
        blaze::Vector<VT3, blaze::columnVector>& xf)
    {
        // Number of integrator steps per simulation step.
        size_t const num_integrator_steps = std::ceil(h / h_max);

        // Actual integrator step
        Real const integrator_step = h / num_integrator_steps;

        ~xf = ~x0;
        for (size_t i = 0; i < num_integrator_steps; ++i)
            (~integrator)(de, t0 + integrator_step * i, integrator_step, ~xf, ~u, ~xf);
    }
}