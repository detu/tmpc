#pragma once

#include <tmpc/integrator/Integrator.hpp>


namespace tmpc
{
    /// @brief CRTP base class for implicit integrators.
    ///
    template <typename Derived>
    class ImplicitIntegrator
    :   public Integrator<Derived>
    {
    protected:
        ImplicitIntegrator() = default;
        ImplicitIntegrator(ImplicitIntegrator const&) = default;
        ImplicitIntegrator(ImplicitIntegrator&&) = default;
        ImplicitIntegrator& operator=(ImplicitIntegrator const&) = default;
        ImplicitIntegrator& operator=(ImplicitIntegrator&&) = default;
        ~ImplicitIntegrator() = default;
    };


    template <
        typename I,
        typename DAE,
        typename DAE_S,
        typename Real,
        typename VT1,
        typename MT1, bool SO1,
        typename VT2,
        typename VT3,
        typename MT2, bool SO2>
    inline void integrate(
        ImplicitIntegrator<I> const& integrator,
        DAE const& dae,
        DAE_S const& dae_s,
        Real t0, Real h, size_t num_integrator_steps, 
        blaze::Vector<VT1, blaze::columnVector> const& x0,
        blaze::Matrix<MT1, SO1> const& Sx,
        blaze::Vector<VT2, blaze::columnVector> const& u,
        blaze::Vector<VT3, blaze::columnVector>& xf,
        blaze::Matrix<MT2, SO2>& Sf)
    {
        // Actual integrator step
        Real const integrator_step = h / num_integrator_steps;

        ~xf = ~x0;
        ~Sf = ~Sx;

        for (size_t i = 0; i < num_integrator_steps; ++i)
            (~integrator)(dae, dae_s, t0 + integrator_step * i, integrator_step, ~xf, ~Sf, ~u, ~xf, ~Sf);
    }


    template <
        typename I,
        typename DAE,
        typename DAE_S,
        typename Residual,
        typename Real,
        typename VT1,
        typename MT1, bool SO1,
        typename VT2,
        typename VT3,
        typename MT2, bool SO2,
        typename VT4,
        typename MT3, bool SO3
    >
    inline void integrate(
        ImplicitIntegrator<I> const& integrator,
        DAE const& dae,
        DAE_S const& dae_s,
        Residual const& res, 
        Real t0, Real h, size_t num_integrator_steps, 
        blaze::Vector<VT1, blaze::columnVector> const& x0,
        blaze::Matrix<MT1, SO1>& S,
        blaze::Vector<VT2, blaze::columnVector> const& u,
        blaze::Vector<VT3, blaze::columnVector>& xf,
        blaze::Matrix<MT2, SO2>& Sf,
        Real& l,
        blaze::Vector<VT4, blaze::columnVector>& g,
        blaze::Matrix<MT3, SO3>& H)
    {
        // Actual integrator step
        Real const integrator_step = h / num_integrator_steps;

        ~xf = ~x0;
        ~Sf = ~S;

        for (size_t i = 0; i < num_integrator_steps; ++i)
            (~integrator)(dae, dae_s, res, t0 + integrator_step * i, integrator_step, ~xf, ~Sf, ~u, ~xf, ~Sf, l, ~g, ~H);
    }
}