#pragma once

#include <tmpc/ocp/OcpKktResidual.hpp>


namespace tmpc
{
    /// @brief Calculate KKT system residual for a given QP and a given solution.
    ///
    template <typename Qp, typename QpSol, typename Real>
    inline void kktResidual(Qp const& qp, QpSol const& sol, OcpKktResidual<Real>& r)
    {
        auto const& graph = qp.graph();

        for (auto v : graph::vertices(graph))
        {
            blaze::DynamicVector<Real, blaze::columnVector> gx = 
                get(qp.Q(), v) * get(sol.x(), v) + trans(get(qp.S(), v)) * get(sol.u(), v) + get(qp.q(), v);
            
            blaze::DynamicVector<Real, blaze::columnVector> gu =
                get(qp.S(), v) * get(sol.x(), v) + get(qp.R(), v) * get(sol.u(), v) + get(qp.r(), v);

            if (auto e = graph.parentEdge(v))
            {
                OcpEdgeDescriptor ed = *e;
                gx -= get(sol.pi(), ed);
            }

            for (auto e : graph::out_edges(v, graph))
            {
                gx += trans(get(qp.A(), e)) * get(sol.pi(), e);
                gu += trans(get(qp.B(), e)) * get(sol.pi(), e);
            }

            put(r.gx(), v, gx);
            put(r.gu(), v, gu);
        }


        for (auto e : graph::edges(graph))
        {
            auto const src = source(e, graph);
            auto const dst = target(e, graph);

            put(r.c(), e, 
                get(qp.A(), e) * get(sol.x(), src)
                + get(qp.B(), e) * get(sol.u(), src)
                + get(qp.b(), e)
                - get(sol.x(), dst));
        }
    }
}