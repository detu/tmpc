#pragma once

#include "detail/BundlePropertyMap.hpp"

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/core/Graph.hpp>
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
            copyProperty(size_map, make_iterator_property_map(size_.begin(), vertexIndex(graph_)), vertices(graph_));

            // Allocate vertex properties of appropriate size
            vertexProperties_.reserve(num_vertices(graph_));
            for (auto const& sz : size_)
                vertexProperties_.emplace_back(sz);
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
            /*
            auto const N = num_vertices(graph_);

            if (N > 0)
            {
                vertexProperties_[N - 1].P_ = 

                for (size_t i = N; i-- > 0; )
                {

                }
            }
            */
            std::vector<boost::default_color_type> color(num_vertices(graph_));

            depth_first_search(graph_,
                RiccatiBackwardVisitor(*this, qp, sol), 
                make_iterator_property_map(color.begin(), get(vertex_index, graph_)), 
                vertex(0, graph_));

            depth_first_search(graph_,
                RiccatiForwardVisitor(*this, qp, sol), 
                make_iterator_property_map(color.begin(), get(vertex_index, graph_)), 
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
        :   public default_dfs_visitor 
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
                    auto const out_e = out_edges(u, g);

                    if (out_e.size() == 1)
                    {
                        auto const e = out_e.front();
                        auto const v = target(e, g);

                        // Alg 1 line 3
                        blaze::DynamicMatrix<Real> PA = trans(get(ws_.P(), v)) * get(qp_.A(), e);
                        blaze::DynamicMatrix<Real> PB = trans(get(ws_.P(), v)) * get(qp_.B(), e);

                        // Alg 1 line 4
                        blaze::DynamicMatrix<Real> BPA = trans(get(qp_.B(), e)) * PA;
                        blaze::DynamicMatrix<Real> BPB = trans(get(qp_.B(), e)) * PB;

                        // Alg 1 line 5
                        blaze::DynamicMatrix<Real> APA = trans(get(qp_.A(), e)) * PA;
                        
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
                        put(ws_.P(), u, get(qp_.Q(), u) + APA - trans(get(ws_.L(), u)) * get(ws_.L(), u));

                        // Alg 1 line 9
                        put(ws_.P(), u, 0.5 * (get(ws_.P(), u) + trans(get(ws_.P(), u))));

                        // Alg 2 line 3
                        put(ws_.l(), u, get(qp_.r(), u) + trans(get(qp_.B(), e)) * (trans(get(ws_.P(), v)) * get(qp_.b(), e) + get(ws_.p(), v)));
                        trsv(get(ws_.Lambda(), u), get(ws_.l(), u), 'L', 'N', 'N');

                        // Alg 2 line 4
                        put(ws_.p(), u, get(qp_.q(), u) + trans(get(qp_.A(), e)) * (trans(get(ws_.P(), v)) * get(qp_.b(), e) + get(ws_.p(), v)) 
                            - trans(get(ws_.L(), u)) * get(ws_.l(), u));
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
        :   public default_dfs_visitor 
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
                }

                // Alg 2 line 8
                put(sol_.u(), u, -(get(ws_.L(), u) * get(sol_.x(), u) + get(ws_.l(), u)));
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
                put(sol_.x(), v, get(qp_.A(), e) * get(sol_.x(), u) + get(qp_.B(), e) * get(sol_.u(), u) + get(qp_.b(), e));

                // Alg 2 line 10
                put(sol_.pi(), e, trans(get(ws_.P(), v)) * get(sol_.x(), v) + get(ws_.p(), v));

                // std::clog << "pi = " << std::endl << get(sol_.pi(), e) << std::endl;
            }


        private:
            ClassicalRiccati& ws_;
            Qp const& qp_;
            QpSol& sol_;
        };


        auto vertexProperties()
        {
            return make_iterator_property_map(vertexProperties_.begin(), vertexIndex(graph_));
        }


        auto P()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::P_, vertexProperties());
        }


        auto Lambda()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::Lambda_, vertexProperties());
        }


        auto L()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::L_, vertexProperties());
        }


        auto p()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::p_, vertexProperties());
        }


        auto l()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::l_, vertexProperties());
        }


        OcpGraph graph_;
        std::vector<OcpSize> size_;
        std::vector<VertexPropertyBundle> vertexProperties_;
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