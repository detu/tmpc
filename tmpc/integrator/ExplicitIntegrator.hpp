#pragma once

#include <tmpc/integrator/Integrator.hpp>


namespace tmpc
{
    /// @brief CRTP base class for explicit integrators.
    ///
    template <typename Derived>
    class ExplicitIntegrator
    :   public Integrator<Derived>
    {
    protected:
        ExplicitIntegrator() = default;
        ExplicitIntegrator(ExplicitIntegrator const&) = default;
        ExplicitIntegrator(ExplicitIntegrator&&) = default;
        ExplicitIntegrator& operator=(ExplicitIntegrator const&) = default;
        ExplicitIntegrator& operator=(ExplicitIntegrator&&) = default;
        ~ExplicitIntegrator() = default;
    };


    template <
        typename I,
        typename ODE,
        typename Real,
        typename VT1,
        typename MT1, bool SO1,
        typename VT2,
        typename VT3,
        typename MT2, bool SO2>
    inline void integrate(
        ExplicitIntegrator<I> const& integrator,
        ODE const& ode,
        Real t0, Real h, Real h_max, 
        blaze::Vector<VT1, blaze::columnVector> const& x0,
        blaze::Matrix<MT1, SO1>& S,
        blaze::Vector<VT2, blaze::columnVector> const& u,
        blaze::Vector<VT3, blaze::columnVector>& xf,
        blaze::Matrix<MT2, SO2>& Sf)
    {
        // Number of integrator steps per simulation step.
        size_t const num_integrator_steps = std::ceil(h / h_max);

        // Actual integrator step
        Real const integrator_step = h / num_integrator_steps;

        ~xf = ~x0;
        ~Sf = ~S;

        for (size_t i = 0; i < num_integrator_steps; ++i)
            (~integrator)(ode, t0 + integrator_step * i, integrator_step, ~xf, ~Sf, ~u, ~xf, ~Sf);
    }


    template <
        typename I,
        typename ODE,
        typename Real,
        typename VT1,
        typename MT1, bool SO1,
        typename VT2,
        typename VT3,
        typename MT2, bool SO2,
        typename VT4,
        typename MT3, bool SO3>
    inline void integrate(
        ExplicitIntegrator<I> const& integrator,
        ODE const& ode,
        Real t0, Real h, Real h_max, 
        blaze::Vector<VT1, blaze::columnVector> const& x0,
        blaze::Matrix<MT1, SO1>& S,
        blaze::Vector<VT2, blaze::columnVector> const& u,
        blaze::Vector<VT3, blaze::columnVector>& xf,
        blaze::Matrix<MT2, SO2>& Sf,
        Real& l,
        blaze::Vector<VT4, blaze::columnVector>& g,
        blaze::Matrix<MT3, SO3>& H)
    {
        // Number of integrator steps per simulation step.
        size_t const num_integrator_steps = std::ceil(h / h_max);

        // Actual integrator step
        Real const integrator_step = h / num_integrator_steps;

        ~xf = ~x0;
        ~Sf = ~S;

        for (size_t i = 0; i < num_integrator_steps; ++i)
            (~integrator)(ode, t0 + integrator_step * i, integrator_step, ~xf, ~Sf, ~u, ~xf, ~Sf, l, ~g, ~H);
    }
}