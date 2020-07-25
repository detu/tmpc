#pragma once

#include <tmpc/ocp/OcpTree.hpp>


namespace tmpc
{
    /// @brief Defines the OCP QP concept.
    ///
    template <typename QP>
    concept OcpQp = requires(QP qp, QP qp1, OcpVertex v, OcpEdge e)
    {
        // Real type
        typename QP::Real;

        // Default constructor
        QP {};

        // Size
        qp.size();

        // QP properties
        qp.graph();
        qp.Q(v);
        qp.R(v);
        qp.S(v);
        qp.q(v);
        qp.r(v);
        qp.C(v);
        qp.D(v);
        qp.ld(v);
        qp.ud(v);
        qp.lx(v);
        qp.ux(v);
        qp.lu(v);
        qp.uu(v);
        qp.A(e);
        qp.B(e);
        qp.b(e);
    };


    /***
     * @brief Assignment of OCP QPs of arbitrary type.
     */
    template <OcpQp Qp1, OcpQp Qp2>
    inline void assign(Qp1& lhs, Qp2 const& rhs)
    {
        if (lhs.size() != rhs.size())
            TMPC_THROW_EXCEPTION(std::invalid_argument("OcpQp of different size cannot be assigned"));

        for (auto v : vertices(lhs.graph()))
        {
            lhs.Q(v, rhs.Q(v));
            lhs.q(v, rhs.q(v));
            lhs.C(v, rhs.C(v));
            lhs.ld(v, rhs.ld(v));
            lhs.ud(v, rhs.ud(v));
            lhs.lx(v, rhs.lx(v));
            lhs.ux(v, rhs.ux(v));
        }

        for (auto v : lhs.graph().branchVertices())
        {
            lhs.R(v, rhs.R(v));
            lhs.S(v, rhs.S(v));
            lhs.r(v, rhs.r(v));
            lhs.D(v, rhs.D(v));
            lhs.lu(v, rhs.lu(v));
            lhs.uu(v, rhs.uu(v));
        }

        for (auto e : edges(lhs.graph()))
        {
            lhs.A(e, rhs.A(e));
            lhs.B(e, rhs.B(e));
            lhs.b(e, rhs.b(e));
        }
    }
}