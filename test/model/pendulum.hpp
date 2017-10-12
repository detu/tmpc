#pragma once

#include <generated_pendulum.h>

#include <tmpc/casadi_interface/GeneratedFunction.hpp>
#include <tmpc/Matrix.hpp>


namespace tmpc :: testing :: model_pendulum
{
    template <typename Kernel>
    class ExplicitOde
    {
    public:
        
        /// \brief Evaluates explicit ODE with sensitivities.
        void operator()(double t, DynamicVector<Kernel> const& x0, DynamicVector<Kernel> const& u0,	DynamicVector<Kernel>& xdot, 
            DynamicMatrix<Kernel>& A, DynamicMatrix<Kernel>& B) const
        {
            ode_({x0.data(), u0.data()}, {xdot.data(), A.data(), B.data()});
        }
    
        /*
        /// \brief Evaluates ODE and quadrature.
        void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
            QuadVector& q, QuadStateMatrix& qA, QuadInputMatrix& qB) const
        {
            _ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data(), q.data(), qA.data(), qB.data(), nullptr, nullptr, nullptr});
        }
    
        /// \brief Evaluates ODE, quadrature and residuals.
        void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
            QuadVector& q, QuadStateMatrix& qA, QuadInputMatrix& qB, ResVector& r, ResStateMatrix& rA, ResInputMatrix& rB) const
        {
            _ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data(), q.data(), qA.data(), qB.data(), r.data(), rA.data(), rB.data()});
        }
    
        void operator()(double t, StateVector const& x0, InputVector const& u0,
                StateVector const& x0_seed, InputVector const& u_seed, StateVector& xdot, StateVector& xdot_sens) const
        {
            static casadi_interface::GeneratedFunction const _ode(pendulum_ode_sens_functions());
            _ode({&t, x0.data(), u0.data(), x0_seed.data(), u_seed.data()}, {xdot.data(), xdot_sens.data()});
        }
    
        /// \brief Evaluates ODE without sensitivities.
        StateVector operator()(double t, StateVector const& x0, InputVector const& u0) const
        {
            StateVector xdot;
            _ode({&t, x0.data(), u0.data()}, {xdot.data(), nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr});
    
            return xdot;
        }
        */

    private:
        casadi_interface::GeneratedFunction const ode_ {pendulum_ode_functions()};
    };
}