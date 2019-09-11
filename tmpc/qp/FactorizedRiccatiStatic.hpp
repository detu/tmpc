#pragma once

#include "detail/BundlePropertyMap.hpp"

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/graph/Graph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/Traits.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/math/Llh.hpp>

#include <blaze/Math.h>

#include <vector>


namespace tmpc
{
    /// @brief Implements factorized Riccati algorithm from [Frison2013]
    /// for static matrix sizes.
    ///
    template <typename Real, size_t NX, size_t NU>
    class FactorizedRiccatiStatic
    {
    public:
        FactorizedRiccatiStatic(OcpGraph const& g)
        :   graph_(g)
        ,   ABPBA_(num_edges(g))
        ,   LL_(num_vertices(g))
        ,   p_(num_vertices(g))
        ,   l_(num_vertices(g))
        ,   LBLA_(num_edges(g))
        ,   Pb_p_(num_edges(g))
        {
        }


        auto size() const
        {
            return make_iterator_property_map(size_.begin(), vertexIndex(graph_));
        }


        auto const& graph() const
        {
            return graph_;
        }


        template <typename Qp, typename QpSol>
        void operator()(Qp const& qp, QpSol& sol)
        {
            for (OcpVertexDescriptor v = num_vertices(graph_); v-- > 0; )
                backwardPassVertex(v, qp);

            for (OcpVertexDescriptor v = 0; v < num_vertices(graph_); ++v)
            {
                if (auto const parent_edge = graph_.parentEdge(v))
                    forwardPassEdge(*parent_edge, qp, sol);

                forwardPassVertex(v, sol);
            }
        }


    private:
        template <typename Qp>
        void backwardPassVertex(OcpVertexDescriptor u, Qp const& qp)
        {
            auto& LL = LL_[u];
            auto Lambda = blaze::submatrix<0, 0, NU, NU>(LL);
            auto L_trans = blaze::submatrix<NU, 0, NX, NU>(LL);   
            auto Lcal = blaze::submatrix<NU, NU, NX, NX>(LL);

            if (out_degree(u, graph_) == 0)
            {
                // Alg 3 line 1
                tmpc::llh(get(qp.Q(), u), Lcal);

                // Alg 2 line 1
                p_[u] = get(qp.q(), u);
            }
            else
            {
                auto const out_e = graph::out_edges(u, graph_);

                if (out_e.size() == 1)
                {
                    auto const e = out_e.front();
                    auto const e_index = get(graph::edge_index, graph_, e);
                    auto const v = target(e, graph_);

                    auto& LBLA = LBLA_[e_index];
                    auto& ABPBA = ABPBA_[e_index];
                    auto& Pb_p = Pb_p_[e_index];
                    auto& l = l_[u];
                    auto& p = p_[u];

                    // Alg 3 line 3
                    auto const Lcal_next = blaze::submatrix<NU, NU, NX, NX>(LL_[v]);
                    blaze::submatrix<0, 0, NX, NU>(LBLA) = trans(Lcal_next) * get(qp.B(), e);
                    blaze::submatrix<0, NU, NX, NX>(LBLA) = trans(Lcal_next) * get(qp.A(), e);

                    // Alg 3 line 4
                    ABPBA = declsym(trans(LBLA) * LBLA);

                    // Alg 3 line 5
                    blaze::submatrix<0, 0, NU, NU>(ABPBA) += declsym(get(qp.R(), u));
                    blaze::submatrix<NU, 0, NX, NU>(ABPBA) += trans(get(qp.S(), u));   
                    blaze::submatrix<NU, NU, NX, NX>(ABPBA) += declsym(get(qp.Q(), u));
                    
                    // llh() or potrf()?
                    // TODO: llh() can be used with adaptors. See if using blaze::SymmetricMatrix improves the performance.
                    // https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cholesky-decomposition
                    tmpc::llh(ABPBA, LL);

                    // Alg 2 line 3.
                    // Pb_p = P_{n+1}^T * b_n + p_{n+1} = \mathcal{L}_{n+1} * \mathcal{L}_{n+1}^T * b_n + p_{n+1}
                    Pb_p = p_[v] + Lcal_next * (trans(Lcal_next) * get(qp.b(), e));
                    l = get(qp.r(), u) + trans(get(qp.B(), e)) * Pb_p;
                    l = inv(Lambda) * l;
                    
                    // Alg 2 line 4
                    p = get(qp.q(), u) + trans(get(qp.A(), e)) * Pb_p - L_trans * l;
                }
                else
                {
                    throw std::invalid_argument("FactorizedRiccatiStatic solver is not implemented on tree QPs yet");
                }
            }
        }


        
        template <typename QpSol>
        void forwardPassVertex(OcpVertexDescriptor u, QpSol& sol) const
        {
            auto const& LL = LL_[u];
            auto const Lambda = blaze::submatrix<0, 0, NU, NU>(LL);
            auto const L_trans = blaze::submatrix<NU, 0, NX, NU>(LL);
            auto const Lcal = blaze::submatrix<NU, NU, NX, NX>(LL);

            if (in_degree(u, graph_) == 0)
            {
                // Root vertex.
                
                // Solve P*x+p=0 by using Cholesky factor of P:
                // \mathcal{L}*(\mathcal{L}^T*x)=-p
                // TODO: this should become faster when the following feature is implemented:
                // https://bitbucket.org/blaze-lib/blaze/issues/284/solving-a-linear-system-with-triangular
                put(sol.x(), u, inv(trans(Lcal)) * (inv(Lcal) * (-p_[u])));
            }

            // Alg 2 line 8
            if (out_degree(u, graph_) > 0)
            {
                // Only non-leaf edges have u.
                //
                // NOTE: need a temporary StaticVector here, otherwise with -O3 we get a SEGFAULT
                // because of improperly aligned vmovapd: 
                // vmovapd 0x8(%rdx),%ymm7
                // 
                // Possibly a compiler bug.
                blaze::StaticVector<Real, NU> const tmp = trans(L_trans) * get(sol.x(), u);
                put(sol.u(), u, -inv(trans(Lambda)) * (l_[u] + tmp));
            }

            /*
            std::clog << "u = " << std::endl << get(sol.u(), u) << std::endl;
            std::clog << "x = " << std::endl << get(sol.x(), u) << std::endl;
            */
        }


        template <typename Qp, typename QpSol>
        void forwardPassEdge(OcpEdgeDescriptor e, Qp const& qp, QpSol& sol) const
        {
            auto const u = source(e, graph_);
            auto const v = target(e, graph_);
            auto const& LL_next = LL_[v];
            auto const& Lcal_next = blaze::submatrix<NU, NU, NX, NX>(LL_next);

            // Alg 2 line 9
            put(sol.x(), v, get(qp.b(), e) + get(qp.A(), e) * get(sol.x(), u) + get(qp.B(), e) * get(sol.u(), u));

            // Alg 2 line 10
            put(sol.pi(), e, p_[v] + Lcal_next * (trans(Lcal_next) * get(sol.x(), v)));

            // std::clog << "pi = " << std::endl << get(sol.pi(), e) << std::endl;
        }


        OcpGraph graph_;
        std::vector<OcpSize> size_;

        // Forcing all matrices to be column major, because of this issue:
        // https://bitbucket.org/blaze-lib/blaze/issues/216
        //
        // Also from [Frison2013]:
        // "Particular attention is given in accessing contiguous data
        // in memory: since all matrices are stored in column-major (or
        // Fortran-like) order, the better performance in matrix-matrix
        // multiplications is obtained when the left matrix is transposed
        // and the right one is not."
        std::vector<blaze::StaticVector<Real, NX>> p_;

        // ABPBA = [B'*P*B, B'*P*A; 
        //        A'*P*B, A'*P*A]
        std::vector<blaze::SymmetricMatrix<blaze::StaticMatrix<Real, NX + NU, NX + NU, blaze::columnMajor>>> ABPBA_;

        // LL = [\Lambda, L;
        //       L', \mathcal{L}]
        std::vector<blaze::LowerMatrix<blaze::StaticMatrix<Real, NX + NU, NX + NU, blaze::columnMajor>>> LL_;
        std::vector<blaze::StaticVector<Real, NU>> l_;

        // LBLA = [L'*B, L'*A]
        std::vector<blaze::StaticMatrix<Real, NX, NU + NX, blaze::columnMajor>> LBLA_;

        std::vector<blaze::StaticVector<Real, NX>> Pb_p_;
    };


    template <typename Real, size_t NX, size_t NU>
    struct RealOf<FactorizedRiccatiStatic<Real, NX, NU>>
    {
        using type = Real;
    };
}