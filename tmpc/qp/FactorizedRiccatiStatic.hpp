#pragma once

#include "detail/BundlePropertyMap.hpp"

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/graph/Graph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/Traits.hpp>
#include <tmpc/Math.hpp>

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
            std::vector<boost::default_color_type> color(num_vertices(graph_));

            depth_first_search(graph_,
                RiccatiBackwardVisitor<Qp, QpSol>(*this, qp, sol), 
                make_iterator_property_map(color.begin(), get(graph::vertex_index, graph_)), 
                vertex(0, graph_));

            depth_first_search(graph_,
                RiccatiForwardVisitor<Qp, QpSol>(*this, qp, sol), 
                make_iterator_property_map(color.begin(), get(graph::vertex_index, graph_)), 
                vertex(0, graph_));
        }


    private:
        template <typename Qp, typename QpSol>
        class RiccatiBackwardVisitor
        :   public graph::default_dfs_visitor 
        {
        public:
            RiccatiBackwardVisitor(FactorizedRiccatiStatic& ws, Qp const& qp, QpSol& sol)
            :   ws_(ws)
            ,   qp_(qp)
            ,   sol_(sol)
            {
            }
        
            
            void finish_vertex(OcpVertexDescriptor u, OcpGraph const& g) const
            {
                auto const& sz_u = get(ws_.size(), u);
                auto& LL = get(ws_.LL(), u);
                auto Lambda = blaze::submatrix<0, 0, NU, NU>(LL);
                auto L_trans = blaze::submatrix<NU, 0, NX, NU>(LL);   
                auto Lcal = blaze::submatrix<NU, NU, NX, NX>(LL);

                if (out_degree(u, g) == 0)
                {
                    // Alg 3 line 1
                    llh(get(qp_.Q(), u), Lcal);

                    // Alg 2 line 1
                    put(ws_.p(), u, get(qp_.q(), u));
                }
                else
                {
                    auto const out_e = graph::out_edges(u, g);

                    if (out_e.size() == 1)
                    {
                        auto const e = out_e.front();
                        auto const v = target(e, g);

                        auto& LBLA = get(ws_.LBLA(), e);
                        auto& ABPBA = get(ws_.ABPBA(), e);
                        auto& Pb_p = get(ws_.Pb_p(), e);
                        auto& l = get(ws_.l(), u);
                        auto& p = get(ws_.p(), u);

                        // Alg 3 line 3
                        auto const Lcal_next = blaze::submatrix<NU, NU, NX, NX>(get(ws_.LL(), v));
                        blaze::submatrix<0, 0, NX, NU>(LBLA) = trans(Lcal_next) * get(qp_.B(), e);
                        blaze::submatrix<0, NU, NX, NX>(LBLA) = trans(Lcal_next) * get(qp_.A(), e);

                        // Alg 3 line 4
                        ABPBA = declsym(trans(LBLA) * LBLA);

                        // Alg 3 line 5
                        blaze::submatrix<0, 0, NU, NU>(ABPBA) += declsym(get(qp_.R(), u));
                        blaze::submatrix<NU, 0, NX, NU>(ABPBA) += trans(get(qp_.S(), u));   
                        blaze::submatrix<NU, NU, NX, NX>(ABPBA) += declsym(get(qp_.Q(), u));
                        
                        // llh() or potrf()?
                        // TODO: llh() can be used with adaptors. See if using blaze::SymmetricMatrix improves the performance.
                        // https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cholesky-decomposition
                        llh(ABPBA, LL);
                        // potrf(LL, 'L');

                        // Alg 2 line 3.
                        // Pb_p = P_{n+1}^T * b_n + p_{n+1} = \mathcal{L}_{n+1} * \mathcal{L}_{n+1}^T * b_n + p_{n+1}
                        Pb_p = get(ws_.p(), v) + Lcal_next * (trans(Lcal_next) * get(qp_.b(), e));
                        l = get(qp_.r(), u) + trans(get(qp_.B(), e)) * Pb_p;
                        l = inv(Lambda) * l;
                        
                        // Alg 2 line 4
                        p = get(qp_.q(), u) + trans(get(qp_.A(), e)) * Pb_p - L_trans * l;
                    }
                    else
                    {
                        throw std::invalid_argument("FactorizedRiccatiStatic solver is not implemented on tree QPs yet");
                    }
                }

                // std::clog << "Lambda = " << std::endl << get(ws_.Lambda(), u) << std::endl;
                // std::clog << "L = " << std::endl << get(ws_.L(), u) << std::endl;
                // std::clog << "l = " << std::endl << get(ws_.l(), u) << std::endl;
                // std::clog << "P = " << std::endl << get(ws_.P(), u) << std::endl;
                // std::clog << "p = " << std::endl << get(ws_.p(), u) << std::endl;
            }


        private:
            FactorizedRiccatiStatic& ws_;
            Qp const& qp_;
            QpSol& sol_;
        };


        template <typename Qp, typename QpSol>
        class RiccatiForwardVisitor
        :   public graph::default_dfs_visitor 
        {
        public:
            RiccatiForwardVisitor(FactorizedRiccatiStatic& ws, Qp const& qp, QpSol& sol)
            :   ws_(ws)
            ,   qp_(qp)
            ,   sol_(sol)
            {
            }
        
            
            void discover_vertex(OcpVertexDescriptor u, OcpGraph const& g) const
            {
                auto const& LL = get(ws_.LL(), u);
                auto const Lambda = blaze::submatrix<0, 0, NU, NU>(LL);
                auto const L_trans = blaze::submatrix<NU, 0, NX, NU>(LL);
                auto const Lcal = blaze::submatrix<NU, NU, NX, NX>(LL);

                if (in_degree(u, g) == 0)
                {
                    // Root vertex.
                    
                    // Solve P*x+p=0 by using Cholesky factor of P:
                    // \mathcal{L}*(\mathcal{L}^T*x)=-p
                    // TODO: this should become faster when the following feature is implemented:
                    // https://bitbucket.org/blaze-lib/blaze/issues/284/solving-a-linear-system-with-triangular
                    put(sol_.x(), u, inv(trans(Lcal)) * (inv(Lcal) * (-get(ws_.p(), u))));
                }

                // Alg 2 line 8
                if (out_degree(u, g) > 0)
                {
                    // Only non-leaf edges have u.
                    put(sol_.u(), u, -inv(trans(Lambda)) * (get(ws_.l(), u) + trans(L_trans) * get(sol_.x(), u)));
                }

                /*
                std::clog << "u = " << std::endl << get(sol_.u(), u) << std::endl;
                std::clog << "x = " << std::endl << get(sol_.x(), u) << std::endl;
                */
            }


            void tree_edge(OcpEdgeDescriptor e, OcpGraph const& g) const
            {
                auto const u = source(e, g);
                auto const v = target(e, g);
                auto const& LL_next = get(ws_.LL(), v);
                auto const& Lcal_next = blaze::submatrix<NU, NU, NX, NX>(LL_next);

                // Alg 2 line 9
                put(sol_.x(), v, get(qp_.b(), e) + get(qp_.A(), e) * get(sol_.x(), u) + get(qp_.B(), e) * get(sol_.u(), u));

                // Alg 2 line 10
                put(sol_.pi(), e, get(ws_.p(), v) + Lcal_next * (trans(Lcal_next) * get(sol_.x(), v)));

                // std::clog << "pi = " << std::endl << get(sol_.pi(), e) << std::endl;
            }


        private:
            FactorizedRiccatiStatic& ws_;
            Qp const& qp_;
            QpSol& sol_;
        };


        auto ABPBA()
        {
            return make_iterator_property_map(ABPBA_.begin(), get(graph::edge_index, graph_));
        }


        auto LL()
        {
            return make_iterator_property_map(LL_.begin(), vertexIndex(graph_));
        }


        auto p()
        {
            return make_iterator_property_map(p_.begin(), vertexIndex(graph_));
        }


        auto l()
        {
            return make_iterator_property_map(l_.begin(), vertexIndex(graph_));
        }


        auto LBLA()
        {
            return make_iterator_property_map(LBLA_.begin(), get(graph::edge_index, graph_));
        }


        auto Pb_p()
        {
            return make_iterator_property_map(Pb_p_.begin(), get(graph::edge_index, graph_));
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