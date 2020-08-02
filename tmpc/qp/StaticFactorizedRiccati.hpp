#pragma once

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/StaticOcpSize.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/math/Llh.hpp>
#include <tmpc/math/Trsv.hpp>
#include <tmpc/Exception.hpp>

#include <blazefeo/math/dense/Syrk.hpp>
#include <blazefeo/math/dense/Potrf.hpp>
#include <blazefeo/math/dense/Trmm.hpp>
#include <blazefeo/math/dense/Trsv.hpp>

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
            for (auto u : vertices(graph_) | views::reverse)
            {
                auto& vd_u = vertexData_[u];
                auto const oe = out_edges(u, graph_);

                if (oe.size() == 0)
                {
                    // Alg 1 line 3
                    auto Lcal = vd_u.Lcal();

                    // TODO: tmpc::llh() segfaults here for matrices of size 9;
                    // using blaze::llh() as a workaround.
                    blaze::llh(qp.Q(u), Lcal);
                    // blazefeo::potrf(qp.Q(u), derestrict(Lcal));
                }
                else
                {
                    // RSQ_tilde is an alias to LL_;
                    // potrf() is performed in-place.
                    // The upper-triangular part of the matrix is untouched.
                    auto& RSQ_tilde = derestrict(vd_u.LL_);

                    // ----------------------
                    // Process first out edge
                    // ----------------------
                    auto e = oe.begin();

                    // Alg 1 line 7
                    blazefeo::trmmRightLower(1., trans(qp.BA(*e)), vertexData_[target(*e, graph_)].Lcal(), transD_);

                    // Alg 1 line 5, 8
                    blazefeo::syrk_ln(1., transD_, 1., qp.H(u), RSQ_tilde);

                    // ---------------------------
                    // Process remaining out edges
                    // ---------------------------
                    while (++e != oe.end())
                    {
                        // Alg 1 line 7
                        blazefeo::trmmRightLower(1., trans(qp.BA(*e)), vertexData_[target(*e, graph_)].Lcal(), transD_);

                        // Alg 1 line 8
                        blazefeo::syrk_ln(1., transD_, 1., RSQ_tilde, RSQ_tilde);
                    }

                    // Alg 1 line 10
                    blazefeo::potrf(RSQ_tilde, RSQ_tilde);
                }
            }
        }


        template <typename Qp>
        void backwardSubstitution(Qp const& qp)
        {
            for (auto u : vertices(graph_) | views::reverse)
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
                    blazefeo::trsvLeftLower(vd_u.Lambda(), r_tilde, vd_u.l_);
                    
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
                    blazefeo::trsvLeftLower(vd_u.Lcal(), -vd_u.p_, vd_u.p_);
                    decltype(auto) x_u = sol.x(u);
                    blazefeo::trsvLeftUpper(trans(vd_u.Lcal()), vd_u.p_, x_u);
                }

                // Alg 3 line 3
                //
                vd_u.l_ += trans(vd_u.L_trans()) * sol.x(u);
                decltype(auto) u_u = sol.u(u);
                blazefeo::trsvLeftUpper(trans(vd_u.Lambda()), -vd_u.l_, u_u);

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

        // Temporary variables for each vertex used by the algorithm.
        std::vector<VertexData> vertexData_;

        blaze::StaticMatrix<Real, NU + NX, NX, blaze::columnMajor> transD_;
    };
}