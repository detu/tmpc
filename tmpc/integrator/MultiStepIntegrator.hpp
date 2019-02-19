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
        MultiStepIntegrator(size_t nx, size_t nu)
        :   nx_(nx)
        ,   nu_(nu)
        ,   x_(nx)
        ,   Ai_(nx, nx)
        ,   Bi_(nx, nu)
        ,   A_(nx, nx)
        ,   B_(nx, nu)
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


        template <typename Integrator, typename ODE, typename VT1, typename VT2, 
            typename VT3, typename MT1, bool SO1, typename MT2, bool SO2>
		void operator()(Integrator const& integrator, ODE const& ode, Real t0, 
            blaze::Vector<VT1, blaze::columnVector> const& x0, blaze::Vector<VT2, blaze::columnVector> const& u, Real h, Real h_max,
            blaze::Vector<VT3, blaze::columnVector>& x1, blaze::Matrix<MT1, SO1>& A, blaze::Matrix<MT2, SO2>& B) const
		{
			// Number of integrator steps per simulation step.
            size_t const num_integrator_steps = ceil(h / h_max);

            // Actual integrator step
            Real const integrator_step = h / num_integrator_steps;

            x_ = x0;
            A_ = blaze::IdentityMatrix<Real>(nx_);

            // Not using ZeroMatrix here because of this:
            // https://bitbucket.org/blaze-lib/blaze/issues/230
            // B_ = blaze::ZeroMatrix<Real>(nx_, nu_);
            B_ = 0.;

            for (size_t i = 0; i < num_integrator_steps; ++i)
            {
                integrator(ode, t0 + integrator_step * i, x_, ~u, integrator_step, x_, Ai_, Bi_);
                A_ = Ai_ * A_;
                B_ = Ai_ * B_ + Bi_;
            }
	
			~x1 = x_;
            ~A = A_;
            ~B = B_;
		}


    private:
        size_t const nx_;
        size_t const nu_;
        mutable blaze::DynamicVector<Real> x_;
        mutable blaze::DynamicMatrix<Real> Ai_;
        mutable blaze::DynamicMatrix<Real> Bi_;
        mutable blaze::DynamicMatrix<Real> A_;
        mutable blaze::DynamicMatrix<Real> B_;
    };
}