#pragma once

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/ocp/OcpSolution.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/Math.hpp>

#include <vector>


namespace tmpc
{
    /// @brief Implements classical Riccati algorithm from [Frison2013]
    /// adapted to tree-structured QPs.
    ///
    template <typename R>
    class DynamicClassicalRiccati
    {
    public:
        using Real = R;


        template <OcpSize Size>
        DynamicClassicalRiccati(Size const& size)
        :   graph_ {size.graph()}
        ,   size_ {size}
        ,   P_(num_vertices(graph_))
        ,   p_(num_vertices(graph_))
        ,   Lambda_(num_vertices(graph_))
        ,   L_(num_vertices(graph_))
        ,   l_(num_vertices(graph_))
        ,   PA_(num_edges(graph_))
        ,   PB_(num_edges(graph_))
        ,   APA_(num_edges(graph_))
        ,   Pb_p_(num_edges(graph_))
        {
            for (auto v : vertices(graph_))
            {
                P_[v].resize(size_.nx(v), size_.nx(v));
                Lambda_[v].resize(size_.nu(v), size_.nu(v));
                L_[v].resize(size_.nu(v), size_.nx(v));
                p_[v].resize(size_.nx(v));
                l_[v].resize(size_.nu(v));
            }

            
            for (auto e : edges(graph_))
            {
                auto const u = source(e, graph_);
                auto const v = target(e, graph_);
                
                PA_[e].resize(size_.nx(v), size_.nx(u));
                PB_[e].resize(size_.nx(v), size_.nu(u));
                APA_[e].resize(size_.nx(u), size_.nx(u));
                Pb_p_[e].resize(size_.nx(v));
            }
        }


        auto const& graph() const noexcept
        {
            return graph_;
        }


        template <OcpQp Qp, typename QpSol>
        void operator()(Qp const& qp, QpSol& sol) const
        {
            for (auto v : graph_.vertices() | std::views::reverse)
                vertexBackward(v, qp);

            for (auto v : graph_.vertices())
            {
                if (auto e = graph_.parentEdge(v))
                    edgeForward(*e, qp, sol);

                vertexForward(v, sol);
            }
        }


    private:
        template <OcpQp Qp>
        void vertexBackward(OcpVertex u, Qp const& qp) const
        {
            P_[u] = declsym(qp.Q(u));
            p_[u] = qp.q(u);

            if (out_degree(u, graph_) > 0)
            {
                Lambda_[u] = declsym(qp.R(u));
                L_[u] = qp.S(u);
                l_[u] = qp.r(u);
                
                for (auto e : out_edges(u, graph_))
                {
                    auto const v = target(e, graph_);

                    auto& PA = PA_[e];
                    auto& PB = PB_[e];
                    auto& APA = APA_[e];
                    auto& Pb_p = Pb_p_[e];

                    // Alg 1 line 3
                    PA = trans(P_[v]) * qp.A(e);
                    PB = trans(P_[v]) * qp.B(e);

                    // Alg 1 line 5
                    APA = trans(qp.A(e)) * PA;
                    APA = 0.5 * (APA + trans(APA));
                    assert(isSymmetric(APA));
                        
                    // Alg 1 line 6
                    Lambda_[u] += declsym(trans(qp.B(e)) * PB);
                    assert(isSymmetric(Lambda_[u]));

                    // Alg 1 line 7
                    L_[u] += trans(qp.B(e)) * PA;

                    // For some reason the code below sometimes results 
                    // in non-positive definite matrix Lambda:
                    //
                    // P += declsym(trans(qp.A(e)) * PA);
                    // P = 0.5 * (P + trans(P));
                    //
                    // This is why we use APA.
                    //
                    // Alg 1 line 8
                    P_[u] += APA;
                    assert(isSymmetric(P_[u]));

                    // Alg 2 line 3
                    Pb_p = trans(P_[v]) * qp.b(e) + p_[v];
                    l_[u] += trans(qp.B(e)) * Pb_p;
                        
                    // Alg 2 line 4
                    p_[u] += trans(qp.A(e)) * Pb_p;
                }

                // Alg 1 line 6
                // llh() or potrf()?
                // TODO: llh() can be used with adaptors. See if using blaze::SymmetricMatrix improves the performance.
                // https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cholesky-decomposition
                //llh(get(this->Lambda(), u), get(this->Lambda(), u));
                assert(isPositiveDefinite(Lambda_[u]));
                potrf(Lambda_[u], 'L');

                // Alg 1 line 7
                trsm(Lambda_[u], L_[u], CblasLeft, CblasLower, 1.);

                // Alg 1 line 8
                P_[u] -= declsym(trans(L_[u]) * L_[u]);
                assert(isPositiveDefinite(P_[u]));

                // Alg 1 line 9
                //put(P(), u, 0.5 * (P_[u] + trans(P_[u])));

                // Alg 2 line 3
                trsv(Lambda_[u], l_[u], 'L', 'N', 'N');

                // Alg 2 line 4
                p_[u] -= trans(L_[u]) * l_[u];
            }

            // std::clog << "Lambda = " << std::endl << Lambda_[u] << std::endl;
            // std::clog << "L = " << std::endl << L_[u] << std::endl;
            // std::clog << "l = " << std::endl << l_[u] << std::endl;
            // std::clog << "P = " << std::endl << P_[u] << std::endl;
            // std::clog << "p = " << std::endl << p_[u] << std::endl;
        }


        template <OcpSolution QpSol>
        void vertexForward(OcpVertex u, QpSol& sol) const
        {
            if (in_degree(u, graph_) == 0)
            {
                // Root vertex
                sol.x(u) = -inv(P_[u]) * p_[u];

                /*
                potrf(get(this->P(), u), 'L');
                put(sol_.x(), u, -get(this->p(), u));
                trsv(get(this->P(), u), get(sol_.x(), u), 'L', 'N', 'N');
                */
            }

            // Alg 2 line 8
            // put(sol_.u(), u, -get(this->L(), u) * get(sol_.x(), u));
            // get(sol_.u(), u) -= get(this->l(), u);
            // trsv(get(this->Lambda(), u), get(sol_.u(), u), 'L', 'T', 'N');

            blaze::DynamicVector<Real> u_tmp = L_[u] * sol.x(u) + l_[u];
            trsv(Lambda_[u], u_tmp, 'L', 'T', 'N');
            sol.u(u) = -u_tmp;

            /*
            std::clog << "u = " << std::endl << get(sol_.u(), u) << std::endl;
            std::clog << "x = " << std::endl << get(sol_.x(), u) << std::endl;
            */
        }


        template <OcpQp Qp, OcpSolution QpSol>
        void edgeForward(OcpEdge e, Qp const& qp, QpSol& sol) const
        {
            auto const u = source(e, graph_);
            auto const v = target(e, graph_);

            // Alg 2 line 9
            sol.x(v) = qp.A(e) * sol.x(u);
            sol.x(v) += qp.B(e) * sol.u(u);
            sol.x(v) += qp.b(e);

            // Alg 2 line 10
            sol.pi(e) = trans(P_[v]) * sol.x(v) + p_[v];

            // std::clog << "pi = " << std::endl << get(sol_.pi(), e) << std::endl;
        }


        OcpTree graph_;
        DynamicOcpSize size_;

        // Forcing all matrices to be column major, because of this issue:
        // https://bitbucket.org/blaze-lib/blaze/issues/216
        //
        // Also from [Frison2013]:
        // "Particular attention is given in accessing contiguous data
        // in memory: since all matrices are stored in column-major (or
        // Fortran-like) order, the better performance in matrix-matrix
        // multiplications is obtained when the left matrix is transposed
        // and the right one is not."
        mutable std::vector<blaze::SymmetricMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>>> P_;
        mutable std::vector<blaze::DynamicVector<Real>> p_;
        mutable std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> Lambda_;
        mutable std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> L_;
        mutable std::vector<blaze::DynamicVector<Real>> l_;

        mutable std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> PA_;
        mutable std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> PB_;
        mutable std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> APA_;
        mutable std::vector<blaze::DynamicVector<Real>> Pb_p_;
    };
}