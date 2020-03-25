#pragma once

#include <tmpc/property_map/BundlePropertyMap.hpp>

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/graph/Graph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/Traits.hpp>
#include <tmpc/Math.hpp>

#include <blaze/Math.h>

#include <vector>


namespace tmpc
{
    /// @brief Implements factorized Riccati algorithm from [Frison2013]
    ///
    template <typename Real>
    class FactorizedRiccati
    {
    public:
        template <typename SizeMap>
        FactorizedRiccati(OcpGraph const& g, SizeMap size_map)
        :   graph_(g)
        ,   size_(num_vertices(g))
        ,   ABPBA_(num_edges(g))
        ,   LL_(num_vertices(g))
        ,   p_(num_vertices(g))
        ,   l_(num_vertices(g))
        ,   LBLA_(num_edges(g))
        ,   Pb_p_(num_edges(g))
        {
            copyProperty(size_map, make_iterator_property_map(size_.begin(), vertexIndex(graph_)), graph::vertices(graph_));
            
            // vertexProperties_.reserve(num_vertices(graph_));
            // for (auto const& sz : size_)
            //     vertexProperties_.emplace_back(sz);

            for (auto v : graph::vertices(graph_))
            {
                auto const& sz = get(size(), v);
                auto v_id = get(graph::vertex_index, graph_, v);

                LL_[v_id].resize(sz.nx() + sz.nu(), sz.nx() + sz.nu());
                p_[v_id].resize(sz.nx());
                l_[v_id].resize(sz.nu());
            }

            
            for (auto e : graph::edges(graph_))
            {
                auto const sz_u = get(size_map, source(e, g));
                auto const sz_v = get(size_map, target(e, g));
                auto edge_id = get(graph::edge_index, graph_, e);
                
                LBLA_[edge_id].resize(sz_v.nx(), sz_u.nu() + sz_u.nx());
                ABPBA_[edge_id].resize(sz_u.nu() + sz_u.nx());
                Pb_p_[edge_id].resize(sz_v.nx());
            }
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
            RiccatiBackwardVisitor(FactorizedRiccati& ws, Qp const& qp, QpSol& sol)
            :   ws_(ws)
            ,   qp_(qp)
            ,   sol_(sol)
            {
            }
        
            
            void finish_vertex(OcpVertexDescriptor u, OcpGraph const& g) const
            {
                auto const& sz_u = get(ws_.size(), u);
                auto& LL = get(ws_.LL(), u);
                auto Lambda = submatrix(LL, 0, 0, sz_u.nu(), sz_u.nu());
                auto L_trans = submatrix(LL, sz_u.nu(), 0, sz_u.nx(), sz_u.nu());   
                auto Lcal = submatrix(LL, sz_u.nu(), sz_u.nu(), sz_u.nx(), sz_u.nx());

                if (out_degree(u, g) == 0)
                {
                    // Alg 3 line 1
                    Lcal = get(qp_.Q(), u);
                    potrf(Lcal, 'L');

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
                        auto const& sz_v = get(ws_.size(), v);

                        auto& LBLA = get(ws_.LBLA(), e);
                        auto& ABPBA = get(ws_.ABPBA(), e);
                        auto& Pb_p = get(ws_.Pb_p(), e);
                        auto& l = get(ws_.l(), u);
                        auto& p = get(ws_.p(), u);

                        // Alg 3 line 3
                        // TODO: dtrmm
                        auto const Lcal_next = submatrix(get(ws_.LL(), v), sz_v.nu(), sz_v.nu(), sz_v.nx(), sz_v.nx());
                        // submatrix(LBLA, 0, 0, sz_v.nx(), sz_u.nu()) = trans(L_next) * get(qp_.B(), e);
                        // submatrix(LBLA, 0, sz_u.nu(), sz_v.nx(), sz_u.nx()) = trans(L_next) * get(qp_.A(), e);
                        submatrix(LBLA, 0, 0, sz_v.nx(), sz_u.nu()) = get(qp_.B(), e);
                        submatrix(LBLA, 0, sz_u.nu(), sz_v.nx(), sz_u.nx()) = get(qp_.A(), e);
                        trmm(LBLA, trans(Lcal_next), CblasLeft, CblasUpper, 1.);

                        // Alg 3 line 4
                        ABPBA = declsym(trans(LBLA) * LBLA);
                        LL = ABPBA;

                        // Alg 3 line 5
                        Lambda += get(qp_.R(), u);
                        L_trans += trans(get(qp_.S(), u));
                        Lcal += get(qp_.Q(), u);
                        
                        // llh() or potrf()?
                        // TODO: llh() can be used with adaptors. See if using blaze::SymmetricMatrix improves the performance.
                        // https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cholesky-decomposition
                        //llh(get(ws_.Lambda(), u), get(ws_.Lambda(), u));
                        potrf(LL, 'L');

                        // Alg 2 line 3.
                        // Pb_p = P_{n+1}^T * b_n + p_{n+1} = \mathcal{L}_{n+1} * \mathcal{L}_{n+1}^T * b_n + p_{n+1}
                        Pb_p = get(qp_.b(), e);
                        trmv(Pb_p, trans(Lcal_next), CblasUpper);
                        trmv(Pb_p, Lcal_next, CblasLower);
                        Pb_p += get(ws_.p(), v);
                        l = get(qp_.r(), u) + trans(get(qp_.B(), e)) * Pb_p;
                        trsv(Lambda, l, 'L', 'N', 'N');
                        
                        // Alg 2 line 4
                        p = get(qp_.q(), u) + trans(get(qp_.A(), e)) * Pb_p - L_trans * l;
                    }
                    else
                    {
                        throw std::invalid_argument("FactorizedRiccati solver is not implemented on tree QPs yet");
                    }
                }

                // std::clog << "Lambda = " << std::endl << get(ws_.Lambda(), u) << std::endl;
                // std::clog << "L = " << std::endl << get(ws_.L(), u) << std::endl;
                // std::clog << "l = " << std::endl << get(ws_.l(), u) << std::endl;
                // std::clog << "P = " << std::endl << get(ws_.P(), u) << std::endl;
                // std::clog << "p = " << std::endl << get(ws_.p(), u) << std::endl;
            }


        private:
            FactorizedRiccati& ws_;
            Qp const& qp_;
            QpSol& sol_;
        };


        template <typename Qp, typename QpSol>
        class RiccatiForwardVisitor
        :   public graph::default_dfs_visitor 
        {
        public:
            RiccatiForwardVisitor(FactorizedRiccati& ws, Qp const& qp, QpSol& sol)
            :   ws_(ws)
            ,   qp_(qp)
            ,   sol_(sol)
            {
            }
        
            
            void discover_vertex(OcpVertexDescriptor u, OcpGraph const& g) const
            {
                auto const& sz_u = get(ws_.size(), u);
                auto const& LL = get(ws_.LL(), u);
                auto const Lambda = submatrix(LL, 0, 0, sz_u.nu(), sz_u.nu());
                auto const L_trans = submatrix(LL, sz_u.nu(), 0, sz_u.nx(), sz_u.nu());
                auto const Lcal = submatrix(LL, sz_u.nu(), sz_u.nu(), sz_u.nx(), sz_u.nx());

                if (in_degree(u, g) == 0)
                {
                    // Root vertex.
                    
                    // Solve P*x+p=0 by using Cholesky factor of P:
                    // \mathcal{L}*(\mathcal{L}^T*x)=-p
                    put(sol_.x(), u, -get(ws_.p(), u));

                    // Solve \mathcal{L}*z=-p
                    trsv(Lcal, get(sol_.x(), u), 'L', 'N', 'N');

                    // Solve \mathcal{L}^T*x=z
                    trsv(Lcal, get(sol_.x(), u), 'L', 'T', 'N');
                }

                // Alg 2 line 8
                // put(sol_.u(), u, -get(ws_.L(), u) * get(sol_.x(), u));
                // get(sol_.u(), u) -= get(ws_.l(), u);
                // trsv(get(ws_.Lambda(), u), get(sol_.u(), u), 'L', 'T', 'N');
                //
                // TODO: avoid using the temporary variable

                blaze::DynamicVector<Real> u_tmp = trans(L_trans) * get(sol_.x(), u) + get(ws_.l(), u);
                trsv(Lambda, u_tmp, 'L', 'T', 'N');
                put(sol_.u(), u, -u_tmp);

                /*
                std::clog << "u = " << std::endl << get(sol_.u(), u) << std::endl;
                std::clog << "x = " << std::endl << get(sol_.x(), u) << std::endl;
                */
            }


            void tree_edge(OcpEdgeDescriptor e, OcpGraph const& g) const
            {
                auto const u = source(e, g);
                auto const v = target(e, g);
                auto const& sz_v = get(ws_.size(), v);
                auto const& LL_next = get(ws_.LL(), v);
                auto const& Lcal_next = submatrix(LL_next, sz_v.nu(), sz_v.nu(), sz_v.nx(), sz_v.nx());

                // Alg 2 line 9
                // TODO: write in 1 line, should not degrade the performance.
                put(sol_.x(), v, get(qp_.A(), e) * get(sol_.x(), u));
                get(sol_.x(), v) += get(qp_.B(), e) * get(sol_.u(), u);
                get(sol_.x(), v) += get(qp_.b(), e);

                // Alg 2 line 10
                // TODO: avoid using temporary pi
                blaze::DynamicVector<Real> pi = get(sol_.x(), v);
                trmv(pi, trans(Lcal_next), CblasUpper);
                trmv(pi, Lcal_next, CblasLower);
                pi += get(ws_.p(), v);
                put(sol_.pi(), e, pi);

                // std::clog << "pi = " << std::endl << get(sol_.pi(), e) << std::endl;
            }


        private:
            FactorizedRiccati& ws_;
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
        std::vector<blaze::DynamicVector<Real>> p_;

        // ABPBA = [B'*P*B, B'*P*A; 
        //        A'*P*B, A'*P*A]
        std::vector<blaze::SymmetricMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>>> ABPBA_;

        // LL = [\Lambda, L;
        //       L', \mathcal{L}]
        std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> LL_;
        std::vector<blaze::DynamicVector<Real>> l_;

        // LBLA = [L'*B, L'*A]
        std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> LBLA_;

        std::vector<blaze::DynamicVector<Real>> Pb_p_;
    };


    template <typename Real>
    struct KernelOf<FactorizedRiccati<Real>>
    {
        using type = BlazeKernel<Real>;
    };


    template <typename Real>
    struct RealOf<FactorizedRiccati<Real>>
    {
        using type = Real;
    };
}