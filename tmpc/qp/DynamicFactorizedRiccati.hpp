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
    /// @brief Implements factorized Riccati algorithm from [Frison2013]
    ///
    template <typename R>
    class DynamicFactorizedRiccati
    {
    public:
        using Real = R;


        template <OcpSize Size>
        DynamicFactorizedRiccati(Size const& size)
        :   graph_ {size.graph()}
        ,   size_ {size}
        ,   p_(num_vertices(graph_))
        ,   RSQ_tilde_(num_vertices(graph_))
        ,   q_tilde_(num_vertices(graph_))
        ,   r_tilde_(num_vertices(graph_))
        ,   LL_(num_vertices(graph_))
        ,   l_(num_vertices(graph_))
        ,   D_(num_edges(graph_))
        ,   d_(num_edges(graph_))
        {
            // vertexProperties_.reserve(num_vertices(graph_));
            // for (auto const& sz : size_)
            //     vertexProperties_.emplace_back(sz);

            for (auto v : vertices(graph_))
            {
                LL_[v].resize(size_.nx(v) + size_.nu(v), size_.nx(v) + size_.nu(v));
                p_[v].resize(size_.nx(v));
                l_[v].resize(size_.nu(v));
                RSQ_tilde_[v].resize(size_.nu(v) + size_.nx(v));
            }

            
            for (auto e : edges(graph_))
            {
                auto const u = source(e, graph_);
                auto const v = target(e, graph_);
                
                D_[e].resize(size_.nx(v), size_.nu(u) + size_.nx(u));
                d_[e].resize(size_.nx(v));
            }
        }


        auto const& graph() const noexcept
        {
            return graph_;
        }


        template <OcpQp Qp, typename QpSol>
        void operator()(Qp const& qp, QpSol& sol)
        {
            backwardFactorization(qp);
            backwardSubstitution(qp);
            forwardSubstitution(qp, sol);
        }


    private:
        template <OcpQp Qp>
        void backwardFactorization(Qp const& qp) const
        {
            for (auto u : vertices(graph_) | views::reverse)
            {
                auto& LL = LL_[u];
                // auto Lambda = submatrix(LL, 0, 0, size_.nu(u), size_.nu(u));
                // auto L_trans = submatrix(LL, size_.nu(u), 0, size_.nx(u), size_.nu(u));   
                auto Lcal = submatrix(LL, size_.nu(u), size_.nu(u), size_.nx(u), size_.nx(u));

                if (out_degree(u, graph_) == 0)
                {
                    // Alg 1 line 3
                    Lcal = qp.Q(u);
                    potrf(Lcal, 'L');
                }
                else
                {
                    // Alg 1 line 5
                    RSQ_tilde_[u] = qp.H(u);
                        
                    for (auto e : out_edges(u, graph_))
                    {
                        auto const v = target(e, graph_);
                        auto& D = D_[e];

                        // Alg 1 line 7
                        auto const Lcal_next = submatrix(LL_[v], size_.nu(v), size_.nu(v), size_.nx(v), size_.nx(v));
                        // submatrix(D, 0, 0, size_.nx(v), size_.nu(u)) = trans(Lcal_next) * qp.B(e);
                        // submatrix(D, 0, size_.nu(u), size_.nx(v), size_.nx(u)) = trans(Lcal_next) * qp.A(e);
                        submatrix(D, 0, 0, size_.nx(v), size_.nu(u)) = qp.B(e);
                        submatrix(D, 0, size_.nu(u), size_.nx(v), size_.nx(u)) = qp.A(e);
                        trmm(D, trans(Lcal_next), CblasLeft, CblasUpper, 1.);

                        // Alg 1 line 8
                        RSQ_tilde_[u] += declsym(trans(D) * D);  
                    }

                    // Alg 1 line 10
                    // llh() or potrf()?
                    // TODO: llh() can be used with adaptors. See if using blaze::SymmetricMatrix improves the performance.
                    // https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cholesky-decomposition
                    blaze::llh(RSQ_tilde_[u], LL_[u]);
                }
            }
        }


        template <OcpQp Qp>
        void backwardSubstitution(Qp const& qp) const
        {
            for (auto u : vertices(graph_) | views::reverse)
            {
                auto& LL = LL_[u];
                auto Lambda = submatrix(LL, 0, 0, size_.nu(u), size_.nu(u));
                auto L_trans = submatrix(LL, size_.nu(u), 0, size_.nx(u), size_.nu(u));

                if (out_degree(u, graph_) == 0)
                {
                    // Alg 2 line 3
                    p_[u] = qp.q(u);
                }
                else
                {
                    // Alg 2 line 5
                    r_tilde_[u] = qp.r(u);
                    q_tilde_[u] = qp.q(u);

                    for (auto e : out_edges(u, graph_))
                    {
                        auto const v = target(e, graph_);
                        auto const Lcal_next = submatrix(LL_[v], size_.nu(v), size_.nu(v), size_.nx(v), size_.nx(v));

                        // Alg 2 line 7.
                        // d = P_{n+1}^T * b_n + p_{n+1} = \mathcal{L}_{n+1} * \mathcal{L}_{n+1}^T * b_n + p_{n+1}
                        d_[e] = qp.b(e);
                        trmv(d_[e], trans(Lcal_next), CblasUpper);
                        trmv(d_[e], Lcal_next, CblasLower);
                        d_[e] += p_[v];

                        // Alg 2 line 8
                        r_tilde_[u] += trans(qp.B(e)) * d_[e];
                        q_tilde_[u] += trans(qp.A(e)) * d_[e];
                }

                    // Alg 2 line 10
                    l_[u] = r_tilde_[u];
                    trsv(Lambda, l_[u], 'L', 'N', 'N');
                    
                    // Alg 2 line 11
                    p_[u] = q_tilde_[u] - L_trans * l_[u];
                }
            }
        }


        template <OcpQp Qp, OcpSolution Solution>
        void forwardSubstitution(Qp const& qp, Solution& sol) const
        {
            for (auto u : graph_.branchVertices())
            {
                auto const& LL = LL_[u];
                auto const Lambda = submatrix(LL, 0, 0, size_.nu(u), size_.nu(u));
                auto const L_trans = submatrix(LL, size_.nu(u), 0, size_.nx(u), size_.nu(u));
                auto const Lcal = submatrix(LL, size_.nu(u), size_.nu(u), size_.nx(u), size_.nx(u));

                if (in_degree(u, graph_) == 0)
                {
                    // Root vertex.
                    
                    // Solve P*x+p=0 by using Cholesky factor of P:
                    // \mathcal{L}*(\mathcal{L}^T*x)=-p
                    decltype(auto) x = sol.x(u);
                    x = -p_[u];
                    
                    // Solve \mathcal{L}*z=-p
                    trsv(Lcal, x, 'L', 'N', 'N');

                    // Solve \mathcal{L}^T*x=z
                    trsv(Lcal, x, 'L', 'T', 'N');
                }

                // Alg 3 line 3
                // put(sol.u(), u, -get(this->L(), u) * get(sol.x(), u));
                // get(sol.u(), u) -= get(this->l(), u);
                // trsv(get(this->Lambda(), u), get(sol.u(), u), 'L', 'T', 'N');
                //
                // TODO: avoid using the temporary variable
                //
                blaze::DynamicVector<Real> u_tmp = trans(L_trans) * sol.x(u) + l_[u];
                trsv(Lambda, u_tmp, 'L', 'T', 'N');
                sol.u(u) = -u_tmp;

                for (auto e : out_edges(u, graph_))
                {
                    auto const v = target(e, graph_);
                    auto const& LL_next = LL_[v];
                    auto const& Lcal_next = submatrix(LL_next, size_.nu(v), size_.nu(v), size_.nx(v), size_.nx(v));

                    // Alg 3 line 5
                    // TODO: write in 1 line, should not degrade the performance.
                    // TODO: avoid using temporary x
                    sol.x(v) = qp.A(e) * sol.x(u);
                    sol.x(v) += qp.B(e) * sol.u(u);
                    sol.x(v) += qp.b(e);

                    // Alg 3 line 6
                    decltype(auto) pi = sol.pi(e);
                    pi = sol.x(v);
                    trmv(pi, trans(Lcal_next), CblasUpper);
                    trmv(pi, Lcal_next, CblasLower);
                    pi += p_[v];
                }
            }
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
        mutable std::vector<blaze::DynamicVector<Real>> p_;

        // RSQ_tilde = [B'*P*B, B'*P*A; 
        //        A'*P*B, A'*P*A]
        mutable std::vector<blaze::SymmetricMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>>> RSQ_tilde_;
        mutable std::vector<blaze::DynamicVector<Real, blaze::columnVector>> q_tilde_;
        mutable std::vector<blaze::DynamicVector<Real, blaze::columnVector>> r_tilde_;

        // LL = [\Lambda, L;
        //       L', \mathcal{L}]
        mutable std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> LL_;
        mutable std::vector<blaze::DynamicVector<Real>> l_;

        // D = [L'*B, L'*A]
        mutable std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> D_;

        mutable std::vector<blaze::DynamicVector<Real>> d_;
    };
}