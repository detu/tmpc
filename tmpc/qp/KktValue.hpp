#pragma once

#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/ocp/OcpSolution.hpp>
#include <tmpc/ocp/OcpKktValue.hpp>


namespace tmpc
{
    /// @brief Calculate KKT system residual for a given QP and a given solution.
    ///
    template <OcpQp Qp, OcpSolution QpSol, OcpKktValue Kkt>
    inline void kktValue(Qp const& qp, QpSol const& sol, Kkt& kkt)
    {
        using Real = typename Kkt::Real;
        
        auto const& graph = qp.graph();

        if (sol.size() != qp.size() || kkt.size() != qp.size())
            TMPC_THROW_EXCEPTION(std::invalid_argument("The qp, sol, and kkt must have the same size"));

        for (auto v : vertices(graph))
        {
            blaze::DynamicVector<Real, blaze::columnVector> gx = 
                qp.Q(v) * sol.x(v) + qp.q(v)
                + sol.lam_ux(v) - sol.lam_lx(v)
                + trans(qp.C(v)) * (sol.lam_ud(v) - sol.lam_ld(v));

            if (out_degree(v, graph) > 0)
                gx += trans(qp.S(v)) * sol.u(v);

            if (auto e = graph.parentEdge(v))
                gx -= sol.pi(*e);

            for (auto e : out_edges(v, graph))
                gx += trans(qp.A(e)) * sol.pi(e);

            kkt.gx(v, gx);
        }

        for (auto v : graph.branchVertices())
        {
            blaze::DynamicVector<Real, blaze::columnVector> gu =
                qp.S(v) * sol.x(v) + qp.R(v) * sol.u(v) + qp.r(v)
                + sol.lam_uu(v) - sol.lam_lu(v)
                + trans(qp.D(v)) * (sol.lam_ud(v) - sol.lam_ld(v));

            for (auto e : out_edges(v, graph))
                gu += trans(qp.B(e)) * sol.pi(e);

            kkt.gu(v, gu);
        }

        for (auto e : edges(graph))
        {
            auto const src = source(e, graph);
            auto const dst = target(e, graph);

            kkt.c(e, 
                qp.A(e) * sol.x(src)
                + qp.B(e) * sol.u(src)
                + qp.b(e)
                - sol.x(dst));
        }
    }
}