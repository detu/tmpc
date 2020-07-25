#pragma once

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/StaticOcpSize.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/math/Llh.hpp>
#include <tmpc/math/SyrkPotrf.hpp>
#include <tmpc/math/Trsv.hpp>
#include <tmpc/Exception.hpp>

#include <vector>


namespace tmpc
{
    /// @brief Implements factorized Riccati algorithm from [Frison2017a]
    /// for static matrix sizes.
    ///
    template <typename Real_, size_t NX_, size_t NU_>
    class StaticFactorizedRiccati
    {
    public:
        using Real = Real_;
        static size_t constexpr NX = NX_;
        static size_t constexpr NU = NU_;


        StaticFactorizedRiccati(OcpTree const& g)
        :   graph_ {g}
        ,   size_ {g}
        ,   vertexData_(num_vertices(g))
        {
        }


        auto const& size() const noexcept
        {
            return size_;
        }


        auto const& graph() const
        {
            return graph_;
        }


        template <typename Qp, typename QpSol>
        void operator()(Qp const& qp, QpSol& sol)
        {
            backwardFactorization(qp);
            backwardSubstitution(qp);
            forwardSubstitution(qp, sol);
        }


    private:
        template <typename Qp>
        void backwardFactorization(Qp const& qp)
        {
            for (auto u : vertices(graph_) | std::views::reverse)
            {
                auto& vd_u = vertexData_[u];

                if (out_degree(u, graph_) == 0)
                {
                    // Alg 1 line 3
                    auto Lcal = vd_u.Lcal();

                    // TODO: tmpc::llh() segfaults here for matrices of size 9;
                    // using blaze::llh() as a workaround.
                    // tmpc::llh(qp.Q(u), Lcal);
                    blaze::llh(qp.Q(u), Lcal);
                }
                else
                {
                    // Alg 1 line 5
                    blaze::SymmetricMatrix<blaze::StaticMatrix<
                        Real, NX + NU, NX + NU, blaze::columnMajor>> RSQ_tilde = qp.H(u);

                    for (auto e : out_edges(u, graph_))
                    {
                        auto const v = target(e, graph_);
                        auto const& vd_v = vertexData_[v];
                        auto const Lcal_next = vd_v.Lcal();

                        // Alg 1 line 7
                        // D = trans(Lcal_next) * get(qp.BA(), e);
                        blaze::StaticMatrix<Real, NX, NU + NX, blaze::columnMajor> D
                            = trans(Lcal_next) * qp.BA(e);

                        // Alg 1 line 8
                        RSQ_tilde += declsym(trans(D) * D);
                    }
                    
                    // Alg 1 line 10
                    blaze::llh(RSQ_tilde, vd_u.LL_);
                }
            }
        }


        template <typename Qp>
        void backwardSubstitution(Qp const& qp)
        {
            for (auto u : vertices(graph_) | std::views::reverse)
            {
                auto& vd_u = vertexData_[u];

                if (out_degree(u, graph_) == 0)
                {
                    // Alg 2 line 3
                    vd_u.p_ = qp.q(u);
                }
                else
                {
                    // Alg 2 line 5
                    blaze::StaticVector<Real, NU> r_tilde = qp.r(u);
                    blaze::StaticVector<Real, NX> q_tilde = qp.q(u);
                                    
                    for (auto e : out_edges(u, graph_))
                    {
                        auto const v = target(e, graph_);
                        auto const& vd_v = vertexData_[v];
                        auto const Lcal_next = vd_v.Lcal();

                        // Alg 2 line 7.
                        // Pb_p = P_v^T * b_e + p_v = \mathcal{L}_v * \mathcal{L}_v^T * b_e + p_v
                        blaze::StaticVector<Real, NX> const tmp1 = trans(Lcal_next) * qp.b(e);
                        blaze::StaticVector<Real, NX> const Pb_p = vd_v.p_ + Lcal_next * tmp1;

                        // Alg 2 line 8.
                        r_tilde += trans(qp.B(e)) * Pb_p;
                        q_tilde += trans(qp.A(e)) * Pb_p;
                    }
                    
                    // Alg 2 line 10
                    vd_u.l_ = inv(vd_u.Lambda()) * r_tilde;
                    
                    // Alg 2 line 11
                    vd_u.p_ = q_tilde - vd_u.L_trans() * vd_u.l_;
                }
            }
        }


        template <typename Qp, typename QpSol>
        void forwardSubstitution(Qp const& qp, QpSol& sol)
        {
            for (auto u : graph_.branchVertices())
            {
                auto& vd_u = vertexData_[u];

                if (in_degree(u, graph_) == 0)
                {
                    // Root vertex.
                    
                    // Solve P*x+p=0 by using Cholesky factor of P:
                    // \mathcal{L}*(\mathcal{L}^T*x)=-p
                    blaze::StaticVector<Real, NX> const tmp1 = inv(vd_u.Lcal()) * vd_u.p_;
                    blaze::StaticVector<Real, NX> const tmp2 = inv(trans(vd_u.Lcal())) * tmp1;

                    sol.x(u) = -tmp2;
                }

                // Alg 3 line 3
                //
                // NOTE: need a temporary StaticVector here, otherwise with -O2 we get a SEGFAULT
                // because of improperly aligned vmovapd: 
                // vmovapd 0x8(%rdx),%ymm7
                // 
                // Possibly a compiler bug.
                // See this issue: https://bitbucket.org/blaze-lib/blaze/issues/290/segfault-in-matrix-matrix-assignment-due
                //
                // NOTE: splitting the x + A * y operation into two and reusing vd_u.l
                // as a temporary, because of the same issue
                //
                blaze::StaticVector<Real, NU>& tmp1 = vd_u.l_;
                tmp1 += trans(vd_u.L_trans()) * sol.x(u);
                blaze::StaticVector<Real, NU> const tmp2 = inv(trans(vd_u.Lambda())) * tmp1;
                sol.u(u) = -tmp2;

                for (auto e : out_edges(u, graph_))
                {
                    auto const v = target(e, graph_);
                    auto const& vd_v = vertexData_[v];

                    // Alg 3 line 5
                    sol.x(v) = qp.b(e) + qp.A(e) * sol.x(u) + qp.B(e) * sol.u(u);

                    // Alg 3 line 6
                    blaze::StaticVector<Real, NX> const tmp = trans(vd_v.Lcal()) * sol.x(v);
                    sol.pi(e) = vd_v.p_ + vd_v.Lcal() * tmp;

                    // std::clog << "pi = " << std::endl << get(sol.pi(), e) << std::endl;
                }
            }
        }


        OcpTree graph_;
        StaticOcpSize<NX, NU, 0> size_;

        struct VertexData
        {
            // Forcing all matrices to be column major, because of this issue:
            // https://bitbucket.org/blaze-lib/blaze/issues/216
            //
            // Also from [Frison2013]:
            // "Particular attention is given in accessing contiguous data
            // in memory: since all matrices are stored in column-major (or
            // Fortran-like) order, the better performance in matrix-matrix
            // multiplications is obtained when the left matrix is transposed
            // and the right one is not."
            
            // LL = [\Lambda, L;
            //       L', \mathcal{L}]
            blaze::LowerMatrix<blaze::StaticMatrix<Real, NX + NU, NX + NU, blaze::columnMajor>> LL_;
            blaze::StaticVector<Real, NU> l_;
            blaze::StaticVector<Real, NX> p_;

            auto Lambda()
            {
                return blaze::submatrix<0, 0, NU, NU>(LL_);
            }


            auto Lambda() const
            {
                return blaze::submatrix<0, 0, NU, NU>(LL_);
            }


            auto L_trans()
            {
                return blaze::submatrix<NU, 0, NX, NU>(LL_);
            }


            auto L_trans() const
            {
                return blaze::submatrix<NU, 0, NX, NU>(LL_);
            }


            auto Lcal()
            {
                return blaze::submatrix<NU, NU, NX, NX>(LL_);
            }


            auto Lcal() const
            {
                return blaze::submatrix<NU, NU, NX, NX>(LL_);
            }
        };

        // Temporary variables for each vertex used by the algorithm
        std::vector<VertexData> vertexData_;
    };
}