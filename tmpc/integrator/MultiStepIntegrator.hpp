#pragma once

#include <blaze/Math.h>


namespace tmpc
{
    /// @brief Integrates on subintervals not bigger than a given length.
    ///
    /// The name of the class might be not the best, because "Multi-Step Integration Methods" are something else,
    /// but I don't have a better name at the moment.
    template <typename Real>
    class MultiStepIntegrator
    {
    public:
        MultiStepIntegrator(size_t nx)
        :   x_(nx)
        {
        }


        template <typename Integrator, typename ODE, typename StateVector0, typename InputVector>
		auto const& operator()(Integrator const& integrator, ODE const& ode, Real t0, StateVector0 const& x0, InputVector const& u, Real h, Real h_max) const
		{
			// Number of integrator steps per simulation step.
            size_t const num_integrator_steps = ceil(h / h_max);

            // Actual integrator step
            Real const integrator_step = h / num_integrator_steps;

            x_ = x0;
            for (size_t i = 0; i < num_integrator_steps; ++i)
                x_ = integrator(ode, t0 + integrator_step * i, x_, u, integrator_step);
	
			return x_;
		}


    private:
        mutable blaze::DynamicVector<Real> x_;
    };
}