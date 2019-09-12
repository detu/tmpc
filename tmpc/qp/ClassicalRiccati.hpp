#pragma once

#include <tmpc/core/detail/BundlePropertyMap.hpp>

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
    /// @brief Implements classical Riccati algorithm from [Frison2013]
    ///
    template <typename Real>
    class ClassicalRiccati
    {
    public:
        template <typename SizeMap>
        ClassicalRiccati(OcpGraph const& g, SizeMap size_map)
        :   graph_(g)
        ,   size_(num_vertices(g))
        ,   P_(num_vertices(g))
        ,   Lambda_(num_vertices(g))
        ,   L_(num_vertices(g))
        ,   p_(num_vertices(g))
        ,   l_(num_vertices(g))
        ,   PA_(num_edges(g))
        ,   PB_(num_edges(g))
        ,   BPA_(num_edges(g))
        ,   BPB_(num_edges(g))
        ,   APA_(num_edges(g))
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

                P_[v_id].resize(sz.nx(), sz.nx());
                Lambda_[v_id].resize(sz.nu(), sz.nu());
                L_[v_id].resize(sz.nu(), sz.nx());
                p_[v_id].resize(sz.nx());
                l_[v_id].resize(sz.nu());
            }

            
            for (auto e : graph::edges(graph_))
            {
                auto const sz_u = get(size_map, source(e, g));
                auto const sz_v = get(size_map, target(e, g));
                auto edge_id = get(graph::edge_index, graph_, e);
                
                PA_[edge_id].resize(sz_v.nx(), sz_u.nx());
                PB_[edge_id].resize(sz_v.nx(), sz_u.nu());
                BPA_[edge_id].resize(sz_u.nu(), sz_u.nx());
                BPB_[edge_id].resize(sz_u.nu(), sz_u.nu());
                APA_[edge_id].resize(sz_u.nx(), sz_u.nx());
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
            RiccatiBackwardVisitor(ClassicalRiccati& ws, Qp const& qp, QpSol& sol)
            :   ws_(ws)
            ,   qp_(qp)
            ,   sol_(sol)
            {
            }
        
            
            void finish_vertex(OcpVertexDescriptor u, OcpGraph const& g) const
            {
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
                        auto& Lambda = get(ws_.Lambda(), u);
                        auto& L = get(ws_.L(), u);
                        auto& P = get(ws_.P(), u);
                        auto& Pb_p = get(ws_.Pb_p(), e);
                        auto& l = get(ws_.l(), u);
                        auto& p = get(ws_.p(), u);

                        // Alg 1 line 3
                        PA = trans(get(ws_.P(), v)) * get(qp_.A(), e);
                        PB = trans(get(ws_.P(), v)) * get(qp_.B(), e);

                        // Alg 1 line 4
                        BPA = trans(get(qp_.B(), e)) * PA;
                        BPB = trans(get(qp_.B(), e)) * PB;
                        assert(isSymmetric(BPB));

                        // Alg 1 line 5
                        APA = trans(get(qp_.A(), e)) * PA;
                        APA = 0.5 * (APA + trans(APA));
                        assert(isSymmetric(APA));
                        
                        // Alg 1 line 6
                        // llh() or potrf()?
                        // TODO: llh() can be used with adaptors. See if using blaze::SymmetricMatrix improves the performance.
                        // https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cholesky-decomposition
                        //llh(get(ws_.Lambda(), u), get(ws_.Lambda(), u));
                        Lambda = declsym(get(qp_.R(), u)) + declsym(BPB);
                        potrf(Lambda, 'L');

                        // Alg 1 line 7
                        L = get(qp_.S(), u) + BPA;

                        if (!isEmpty(get(qp_.S(), u)))  // prevent trsm() from complaining about lda=0
                            trsm(Lambda, L, CblasLeft, CblasLower, 1.);
                        
                        // Alg 1 line 8
                        P = declsym(get(qp_.Q(), u)) + declsym(APA) - declsym(trans(L) * L);

                        // Alg 1 line 9
                        //put(ws_.P(), u, 0.5 * (get(ws_.P(), u) + trans(get(ws_.P(), u))));

                        // Alg 2 line 3
                        Pb_p = trans(get(ws_.P(), v)) * get(qp_.b(), e) + get(ws_.p(), v);
                        l = get(qp_.r(), u) + trans(get(qp_.B(), e)) * Pb_p;
                        trsv(Lambda, l, 'L', 'N', 'N');
                        
                        // Alg 2 line 4
                        p = get(qp_.q(), u) + trans(get(qp_.A(), e)) * Pb_p - trans(L) * l;
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
                // put(sol_.u(), u, -get(ws_.L(), u) * get(sol_.x(), u));
                // get(sol_.u(), u) -= get(ws_.l(), u);
                // trsv(get(ws_.Lambda(), u), get(sol_.u(), u), 'L', 'T', 'N');

                blaze::DynamicVector<Real> u_tmp = get(ws_.L(), u) * get(sol_.x(), u) + get(ws_.l(), u);
                trsv(get(ws_.Lambda(), u), u_tmp, 'L', 'T', 'N');
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

        // Forcing all matrices to be column major, because of this issue:
        // https://bitbucket.org/blaze-lib/blaze/issues/216
        //
        // Also from [Frison2013]:
        // "Particular attention is given in accessing contiguous data
        // in memory: since all matrices are stored in column-major (or
        // Fortran-like) order, the better performance in matrix-matrix
        // multiplications is obtained when the left matrix is transposed
        // and the right one is not."
        std::vector<blaze::SymmetricMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>>> P_;
        std::vector<blaze::DynamicVector<Real>> p_;
        std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> Lambda_;
        std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> L_;
        std::vector<blaze::DynamicVector<Real>> l_;

        std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> PA_;
        std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> PB_;
        std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> BPA_;
        std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> BPB_;
        std::vector<blaze::DynamicMatrix<Real, blaze::columnMajor>> APA_;
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