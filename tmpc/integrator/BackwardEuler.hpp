#pragma once

#include <tmpc/integrator/ButcherTableau.hpp>


namespace tmpc
{
    /// @brief Butcher tableau for backward Euler mathod.
    template <typename Real>
    inline ButcherTableau<Real> backwardEuler()
    {
        using Tableau = ButcherTableau<Real>;
        return Tableau(
            typename Tableau::A_type {{Real(1)}}, 
            typename Tableau::b_type {Real(1)}, 
            typename Tableau::c_type {Real(1)}
        );
    }
}