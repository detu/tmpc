#pragma once

#include <tmpc/ocp/OcpTree.hpp>


namespace tmpc
{
    /// @brief Defines the OCP solution concept.
    ///
    template <typename Solution>
    concept OcpSolution = requires(Solution sol, OcpVertex v, OcpEdge e)
    {
        typename Solution::Real;

        sol.graph();
        sol.size();

        sol.x(v);
        sol.u(v);
        sol.lam_lx(v);
        sol.lam_ux(v);
        sol.lam_lu(v);
        sol.lam_uu(v);
        sol.lam_ld(v);
        sol.lam_ud(v);
        sol.pi(e);
    };
}