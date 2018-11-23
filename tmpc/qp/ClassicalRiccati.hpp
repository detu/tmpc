#pragma once

#include "detail/BundlePropertyMap.hpp"

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/graph/Graph.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/Traits.hpp>

#include <blaze/Math.h>

#include <vector>


namespace tmpc
{
    template <typename Real>
    class ClassicalRiccati
    {
    public:
        template <typename SizeMap>
        ClassicalRiccati(OcpGraph const& g, SizeMap size_map)
        :   graph_(g)
        ,   size_(num_vertices(g))
        {
            copyProperty(size_map, make_iterator_property_map(size_.begin(), vertexIndex(graph_)), graph::vertices(graph_));

            // Allocate vertex properties of appropriate size
            auto const nv = num_vertices(graph_);
            P_.resize(nv);
            Lambda_.resize(nv);
            L_.resize(nv);
            p_.resize(nv);
            l_.resize(nv);
            
            // vertexProperties_.reserve(num_vertices(graph_));
            // for (auto const& sz : size_)
            //     vertexProperties_.emplace_back(sz);

            auto const vertex_id = get(graph::vertex_index, graph_);
            for (auto v : graph::vertices(graph_))
            {
                auto const& sz = get(size(), v);
                auto v_id = get(vertex_id, v);

                P_[v_id].resize(sz.nx(), sz.nx());
                Lambda_[v_id].resize(sz.nu(), sz.nu());
                L_[v_id].resize(sz.nu(), sz.nx());
                p_[v_id].resize(sz.nx());
                l_[v_id].resize(sz.nu());
            }

            // Allocate edge properties of appropriate size
            auto const ne = num_edges(graph_);
            PA_.resize(ne);
            PB_.resize(ne);
            BPA_.resize(ne);
            BPB_.resize(ne);
            APA_.resize(ne);
            Pb_p_.resize(ne);
            
            auto const edge_id = get(graph::edge_index, graph_);
            for (auto e : graph::edges(graph_))
            {
                auto const sz_u = get(size_map, source(e, g));
                auto const sz_v = get(size_map, target(e, g));
                
                PA_[get(edge_id, e)].resize(sz_v.nx(), sz_u.nx());
                PB_[get(edge_id, e)].resize(sz_v.nx(), sz_u.nu());
                BPA_[get(edge_id, e)].resize(sz_u.nu(), sz_u.nx());
                BPB_[get(edge_id, e)].resize(sz_u.nu(), sz_u.nu());
                APA_[get(edge_id, e)].resize(sz_u.nx(), sz_u.nx());
                Pb_p_[get(edge_id, e)].resize(sz_v.nx());
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
                RiccatiBackwardVisitor(*this, qp, sol), 
                make_iterator_property_map(color.begin(), get(graph::vertex_index, graph_)), 
                vertex(0, graph_));

            depth_first_search(graph_,
                RiccatiForwardVisitor(*this, qp, sol), 
                make_iterator_property_map(color.begin(), get(graph::vertex_index, graph_)), 
                vertex(0, graph_));
        }


    private:
        struct VertexPropertyBundle
        {
            VertexPropertyBundle(OcpSize const& sz)
            :   P_(sz.nx(), sz.nx())
            ,   Lambda_(sz.nu(), sz.nu())
            ,   L_(sz.nu(), sz.nx())
            ,   p_(sz.nx())
            ,   l_(sz.nu())
            {
            }

            
            blaze::DynamicMatrix<Real> P_;
            blaze::DynamicMatrix<Real> Lambda_;
            blaze::DynamicMatrix<Real> L_;
            blaze::DynamicVector<Real> p_;
            blaze::DynamicVector<Real> l_;
        };


        template <typename Qp, typename QpSol>
        class RiccatiBackwardVisitor
        :   public graph::default_dfs_visitor 
        {
        public:
            RiccatiBackwardVisitor(ClassicalRiccati& ws, Qp const& qp, QpSol& sol)
            :   ws_(ws)
            ,   qp_(qp)
            ,   sol_(sol)
            {
            }
        
            
            void finish_vertex(OcpVertexDescriptor u, OcpGraph const& g) const
            {
                // std::clog << "Backward, Vertex " << u << std::endl;

                if (out_degree(u, g) == 0)
                {
                    // Alg 1 line 1
                    put(ws_.P(), u, get(qp_.Q(), u));

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

                        auto& PA = get(ws_.PA(), e);
                        auto& PB = get(ws_.PB(), e);
                        auto& BPA = get(ws_.BPA(), e);
                        auto& BPB = get(ws_.BPB(), e);
                        auto& APA = get(ws_.APA(), e);                        

                        // Alg 1 line 3
                        PA = trans(get(ws_.P(), v)) * get(qp_.A(), e);
                        PB = trans(get(ws_.P(), v)) * get(qp_.B(), e);

                        // Alg 1 line 4
                        BPA = trans(get(qp_.B(), e)) * PA;
                        BPB = trans(get(qp_.B(), e)) * PB;

                        // Alg 1 line 5
                        APA = trans(get(qp_.A(), e)) * PA;
                        
                        // Alg 1 line 6
                        // llh() or potrf()?
                        // TODO: llh() can be used with adaptors. See if using blaze::SymmetricMatrix improves the performance.
                        // https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cholesky-decomposition
                        //llh(get(ws_.Lambda(), u), get(ws_.Lambda(), u));
                        put(ws_.Lambda(), u, get(qp_.R(), u) + BPB);
                        potrf(get(ws_.Lambda(), u), 'L');

                        // Alg 1 line 7
                        put(ws_.L(), u, get(qp_.S(), u) + BPA);

                        if (!isEmpty(get(qp_.S(), u)))  // prevent trsm() from complaining about lda=0
                            trsm(get(ws_.Lambda(), u), get(ws_.L(), u), CblasLeft, CblasLower, 1.);
                        
                        // Alg 1 line 8
                        put(ws_.P(), u, get(qp_.Q(), u) + APA);
                        get(ws_.P(), u) -= trans(get(ws_.L(), u)) * get(ws_.L(), u);

                        // Alg 1 line 9
                        //put(ws_.P(), u, 0.5 * (get(ws_.P(), u) + trans(get(ws_.P(), u))));

                        // Alg 2 line 3
                        put(ws_.Pb_p(), e, trans(get(ws_.P(), v)) * get(qp_.b(), e));
                        get(ws_.Pb_p(), e) += get(ws_.p(), v);

                        put(ws_.l(), u, get(qp_.r(), u));
                        get(ws_.l(), u) += trans(get(qp_.B(), e)) * get(ws_.Pb_p(), e);
                        
                        trsv(get(ws_.Lambda(), u), get(ws_.l(), u), 'L', 'N', 'N');

                        // Alg 2 line 4
                        put(ws_.p(), u, get(qp_.q(), u));
                        get(ws_.p(), u) += trans(get(qp_.A(), e)) * get(ws_.Pb_p(), e);
                        get(ws_.p(), u) -= trans(get(ws_.L(), u)) * get(ws_.l(), u);
                    }
                    else
                    {
                        throw std::invalid_argument("ClassicalRiccati solver is not implemented on tree QPs yet");
                    }
                }

                // std::clog << "Lambda = " << std::endl << get(ws_.Lambda(), u) << std::endl;
                // std::clog << "L = " << std::endl << get(ws_.L(), u) << std::endl;
                // std::clog << "l = " << std::endl << get(ws_.l(), u) << std::endl;
                // std::clog << "P = " << std::endl << get(ws_.P(), u) << std::endl;
                // std::clog << "p = " << std::endl << get(ws_.p(), u) << std::endl;
            }


        private:
            ClassicalRiccati& ws_;
            Qp const& qp_;
            QpSol& sol_;
        };


        template <typename Qp, typename QpSol>
        class RiccatiForwardVisitor
        :   public graph::default_dfs_visitor 
        {
        public:
            RiccatiForwardVisitor(ClassicalRiccati& ws, Qp const& qp, QpSol& sol)
            :   ws_(ws)
            ,   qp_(qp)
            ,   sol_(sol)
            {
            }
        
            
            void discover_vertex(OcpVertexDescriptor u, OcpGraph const& g) const
            {
                // std::clog << "Forward, Vertex " << u << std::endl;

                if (in_degree(u, g) == 0)
                {
                    // Root vertex
                    put(sol_.x(), u, -inv(get(ws_.P(), u)) * get(ws_.p(), u));

                    /*
                    potrf(get(ws_.P(), u), 'L');
                    put(sol_.x(), u, -get(ws_.p(), u));
                    trsv(get(ws_.P(), u), get(sol_.x(), u), 'L', 'N', 'N');
                    */
                }

                // Alg 2 line 8
                put(sol_.u(), u, -get(ws_.L(), u) * get(sol_.x(), u));
                get(sol_.u(), u) -= get(ws_.l(), u);
                trsv(get(ws_.Lambda(), u), get(sol_.u(), u), 'L', 'T', 'N');

                /*
                std::clog << "u = " << std::endl << get(sol_.u(), u) << std::endl;
                std::clog << "x = " << std::endl << get(sol_.x(), u) << std::endl;
                */
            }


            void tree_edge(OcpEdgeDescriptor e, OcpGraph const& g) const
            {
                // std::clog << "Forward, Edge " << e << std::endl;
                auto const u = source(e, g);
                auto const v = target(e, g);

                // Alg 2 line 9
                put(sol_.x(), v, get(qp_.A(), e) * get(sol_.x(), u));
                get(sol_.x(), v) += get(qp_.B(), e) * get(sol_.u(), u);
                get(sol_.x(), v) += get(qp_.b(), e);

                // Alg 2 line 10
                put(sol_.pi(), e, trans(get(ws_.P(), v)) * get(sol_.x(), v) + get(ws_.p(), v));

                // std::clog << "pi = " << std::endl << get(sol_.pi(), e) << std::endl;
            }


        private:
            ClassicalRiccati& ws_;
            Qp const& qp_;
            QpSol& sol_;
        };


        template <typename Qp>
        void vertexBackward(OcpVertexDescriptor u, Qp const& qp)
        {
            // std::clog << "Backward, Vertex " << u << std::endl;
            auto const& g = graph_;

            if (out_degree(u, g) == 0)
            {
                // Alg 1 line 1
                P_[u] = get(qp.Q(), u);

                // Alg 2 line 1
                p_[u] = get(qp.q(), u);
            }
            else
            {
                auto const out_e = graph::out_edges(u, g);

                if (out_e.size() == 1)
                {
                    auto const e = out_e.front();
                    auto const e_id = get(get(graph::edge_index, g), e);
                    auto const v = target(e, g);

                    // Alg 1 line 3
                    PA_[e_id] = trans(P_[v]) * get(qp.A(), e);
                    PB_[e_id] = trans(P_[v]) * get(qp.B(), e);

                    // Alg 1 line 4
                    BPA_[e_id] = trans(get(qp.B(), e)) * PA_[e_id];
                    BPB_[e_id] = trans(get(qp.B(), e)) * PB_[e_id];

                    // Alg 1 line 5
                    APA_[e_id] = trans(get(qp.A(), e)) * PA_[e_id];
                    
                    // Alg 1 line 6
                    // llh() or potrf()?
                    // TODO: llh() can be used with adaptors. See if using blaze::SymmetricMatrix improves the performance.
                    // https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cholesky-decomposition
                    //llh(get(Lambda(), u), get(Lambda(), u));
                    Lambda_[u] = get(qp.R(), u) + BPB_[e_id];
                    potrf(Lambda_[u], 'L');

                    // Alg 1 line 7
                    L_[u] = get(qp.S(), u) + BPA_[e_id];

                    if (!isEmpty(L_[u]))  // prevent trsm() from complaining about lda=0
                        trsm(Lambda_[u], L_[u], CblasLeft, CblasLower, 1.);
                    
                    // Alg 1 line 8
                    P_[u] = get(qp.Q(), u) + APA_[e_id];
                    P_[u] -= trans(L_[u]) * L_[u];

                    // Alg 1 line 9
                    P_[u] = 0.5 * (P_[u] + trans(P_[u]));
                    //P_[u] += trans(P_[u]);
                    //P_[u] *= 0.5;

                    // Alg 2 line 3
                    Pb_p_[e_id] = trans(P_[v]) * get(qp.b(), e);
                    Pb_p_[e_id] += p_[v];

                    l_[u] = get(qp.r(), u);
                    l_[u] += trans(get(qp.B(), e)) * Pb_p_[e_id];
                    
                    trsv(Lambda_[u], l_[u], 'L', 'N', 'N');

                    // Alg 2 line 4
                    p_[u] = get(qp.q(), u);
                    p_[u] += trans(get(qp.A(), e)) * Pb_p_[e_id];
                    p_[u] -= trans(L_[u]) * l_[u];
                }
                else
                {
                    throw std::invalid_argument("ClassicalRiccati solver is not implemented on tree QPs yet");
                }
            }

            // std::clog << "Lambda = " << std::endl << get(Lambda(), u) << std::endl;
            // std::clog << "L = " << std::endl << get(L(), u) << std::endl;
            // std::clog << "l = " << std::endl << get(l(), u) << std::endl;
            // std::clog << "P = " << std::endl << get(P(), u) << std::endl;
            // std::clog << "p = " << std::endl << get(p(), u) << std::endl;
        }


        // auto vertexProperties()
        // {
        //     return make_iterator_property_map(vertexProperties_.begin(), vertexIndex(graph_));
        // }


        auto P()
        {
            return make_iterator_property_map(P_.begin(), vertexIndex(graph_));
        }


        auto Lambda()
        {
            return make_iterator_property_map(Lambda_.begin(), vertexIndex(graph_));
        }


        auto L()
        {
            return make_iterator_property_map(L_.begin(), vertexIndex(graph_));
        }


        auto p()
        {
            return make_iterator_property_map(p_.begin(), vertexIndex(graph_));
        }


        auto l()
        {
            return make_iterator_property_map(l_.begin(), vertexIndex(graph_));
        }


        auto PA()
        {
            return make_iterator_property_map(PA_.begin(), get(graph::edge_index, graph_));
        }


        auto PB()
        {
            return make_iterator_property_map(PB_.begin(), get(graph::edge_index, graph_));
        }


        auto BPA()
        {
            return make_iterator_property_map(BPA_.begin(), get(graph::edge_index, graph_));
        }


        auto BPB()
        {
            return make_iterator_property_map(BPB_.begin(), get(graph::edge_index, graph_));
        }


        auto APA()
        {
            return make_iterator_property_map(APA_.begin(), get(graph::edge_index, graph_));
        }


        auto Pb_p()
        {
            return make_iterator_property_map(Pb_p_.begin(), get(graph::edge_index, graph_));
        }


        OcpGraph graph_;
        std::vector<OcpSize> size_;
        //std::vector<VertexPropertyBundle> vertexProperties_;
        std::vector<blaze::DynamicMatrix<Real>> P_;
        std::vector<blaze::DynamicVector<Real>> p_;
        std::vector<blaze::DynamicMatrix<Real>> Lambda_;
        std::vector<blaze::DynamicMatrix<Real>> L_;
        std::vector<blaze::DynamicVector<Real>> l_;

        std::vector<blaze::DynamicMatrix<Real>> PA_;
        std::vector<blaze::DynamicMatrix<Real>> PB_;
        std::vector<blaze::DynamicMatrix<Real>> BPA_;
        std::vector<blaze::DynamicMatrix<Real>> BPB_;
        std::vector<blaze::DynamicMatrix<Real>> APA_;
        std::vector<blaze::DynamicVector<Real>> Pb_p_;
    };


    template <typename Real>
    struct KernelOf<ClassicalRiccati<Real>>
    {
        using type = BlazeKernel<Real>;
    };


    template <typename Real>
    struct RealOf<ClassicalRiccati<Real>>
    {
        using type = Real;
    };
}