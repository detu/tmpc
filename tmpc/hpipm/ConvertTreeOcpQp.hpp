#pragma once

#include <tmpc/hpipm/TreeOcpQp.hpp>
#include <tmpc/qp/OcpQp.hpp>


namespace tmpc :: hpipm
{
    /// @brief Copy QP to HPIPM data structures.
    template <tmpc::OcpQp Qp, typename Real>
    inline void convertQp(Qp const& src, TreeOcpQp<Real>& dst)
    {
        auto const& graph = src.graph();

    #if 1
        for (auto v : vertices(graph))
        {
            dst.set_Q(v, src.Q(v));
            dst.set_q(v, src.q(v));
            dst.set_lbx(v, src.lx(v));
            dst.set_ubx(v, src.ux(v));
            dst.set_C(v, src.C(v));
            dst.set_lg(v, src.ld(v));
            dst.set_ug(v, src.ud(v));

            // **** No soft constraints support yet!
            //
            // dst->set_Zl(v, Zl_[v].data());
            // dst->set_Zu(v, Zu_[v].data());
            // dst->set_zl(v, zl_[v].data());
            // dst->set_zu(v, zu_[v].data());
            // dst->set_idxs(v, idxs_[v].get());
            // dst->set_lls(v, lbls_[v].data());
            // dst->set_lus(v, lbus_[v].data());
        }

        for (auto v : graph.branchVertices())
        {
            // Converting matrices to unpadded format before passing them to HPIPM functions
            dst.set_R(v, src.R(v));
            dst.set_S(v, src.S(v));
            dst.set_r(v, src.r(v));
            dst.set_lbu(v, src.lu(v));
            dst.set_ubu(v, src.uu(v));
            dst.set_D(v, src.D(v));
        }

        for (auto e : edges(graph))
        {
            dst.set_A(e, src.A(e));
            dst.set_B(e, src.B(e));
            dst.set_b(e, src.b(e));
        }
    #else
        using UMatrix = UnpaddedMatrix<Real, blaze::columnMajor>;

        std::vector<UMatrix>
        cvt_colmaj_to_tree_ocp_qp(
            A, B, b,
            Q, S, R, q, r,
            idxb, d_lb, d_ub,
            C, D, d_lg, d_ug,
            Zl, Zu, zl, zu,
            idxs, d_ls, d_us, dst);
    #endif
    }
}