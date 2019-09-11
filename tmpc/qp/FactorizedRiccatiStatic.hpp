#pragma once

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/graph/Graph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/Traits.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/math/Llh.hpp>

#include <boost/throw_exception.hpp>

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
        ,   vertexData_(num_vertices(g))
        {
        }


        auto size() const
        {
            BOOST_THROW_EXCEPTION(std::logic_error("Function not implemented"));

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
            auto& vd_u = vertexData_[u];

            if (out_degree(u, graph_) == 0)
            {
                // Alg 3 line 1
                auto Lcal = vd_u.Lcal();
                tmpc::llh(get(qp.Q(), u), Lcal);

                // Alg 2 line 1
                vd_u.p_ = get(qp.q(), u);
            }
            else
            {
                auto const out_e = graph::out_edges(u, graph_);

                if (out_e.size() == 1)
                {
                    auto const e = out_e.front();
                    auto const v = target(e, graph_);
                    auto const& vd_v = vertexData_[v];

                    // Alg 3 line 3
                    auto const Lcal_next = vd_v.Lcal();
                    // LBLA_ = trans(Lcal_next) * get(qp.BA(), e);
                    blaze::submatrix<0, 0, NX, NU>(LBLA_) = trans(Lcal_next) * get(qp.B(), e);
                    blaze::submatrix<0, NU, NX, NX>(LBLA_) = trans(Lcal_next) * get(qp.A(), e);

                    // Alg 3 line 4, 5
                    {
                        decltype(auto) ABPBA = derestrict(vd_u.LL_);
                        ABPBA = declsym(get(qp.H(), u)) + declsym(trans(LBLA_) * LBLA_);
                        tmpc::llh(ABPBA);
                    }

                    // Alg 2 line 3.
                    // Pb_p = P_{n+1}^T * b_n + p_{n+1} = \mathcal{L}_{n+1} * \mathcal{L}_{n+1}^T * b_n + p_{n+1}
                    blaze::StaticVector<Real, NX> const tmp1 = trans(Lcal_next) * get(qp.b(), e);
                    blaze::StaticVector<Real, NX> const Pb_p = vd_v.p_ + Lcal_next * tmp1;
                        
                    auto& l = vd_u.l_;
                    blaze::StaticVector<Real, NU> const tmp2 = get(qp.r(), u) + trans(get(qp.B(), e)) * Pb_p;
                    l = inv(vd_u.Lambda()) * tmp2;
                    
                    // Alg 2 line 4
                    auto& p = vd_u.p_;
                    p = get(qp.q(), u) + trans(get(qp.A(), e)) * Pb_p - vd_u.L_trans() * l;
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
            auto const& vd_u = vertexData_[u];

            if (in_degree(u, graph_) == 0)
            {
                // Root vertex.
                
                // Solve P*x+p=0 by using Cholesky factor of P:
                // \mathcal{L}*(\mathcal{L}^T*x)=-p
                // TODO: this should become faster when the following feature is implemented:
                // https://bitbucket.org/blaze-lib/blaze/issues/284/solving-a-linear-system-with-triangular
                blaze::StaticVector<Real, NX> const tmp1 = inv(vd_u.Lcal()) * vd_u.p_;
                blaze::StaticVector<Real, NX> const tmp2 = inv(trans(vd_u.Lcal())) * tmp1;
                put(sol.x(), u, -tmp2);
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
                blaze::StaticVector<Real, NU> const tmp1 = vd_u.l_ + trans(vd_u.L_trans()) * get(sol.x(), u);
                blaze::StaticVector<Real, NU> const tmp2 = inv(trans(vd_u.Lambda())) * tmp1;
                put(sol.u(), u, -tmp2);
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
            auto const& vd_v = vertexData_[v];

            // Alg 2 line 9
            put(sol.x(), v, get(qp.b(), e) + get(qp.A(), e) * get(sol.x(), u) + get(qp.B(), e) * get(sol.u(), u));

            // Alg 2 line 10
            blaze::StaticVector<Real, NX> const tmp = trans(vd_v.Lcal()) * get(sol.x(), v);
            put(sol.pi(), e, vd_v.p_ + vd_v.Lcal() * tmp);

            // std::clog << "pi = " << std::endl << get(sol.pi(), e) << std::endl;
        }


        OcpGraph graph_;
        std::vector<OcpSize> size_;

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

        // LBLA = [L'*B, L'*A]
        blaze::StaticMatrix<Real, NX, NU + NX, blaze::columnMajor> LBLA_;
    };


    template <typename Real, size_t NX, size_t NU>
    struct RealOf<FactorizedRiccatiStatic<Real, NX, NU>>
    {
        using type = Real;
    };
}