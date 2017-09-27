#include "pendulum.hpp"

#include <tmpc/EigenKernel.hpp>

namespace tmpc :: testing :: model_pendulum
{
    void f()
    {
        using Kernel = EigenKernel<double>;

        ExplicitOde<Kernel> ode;
        DynamicVector<Kernel> x0, xdot, u0;
        DynamicMatrix<Kernel> A, B;

        ode(0., x0, u0, xdot, A, B);
    }
}