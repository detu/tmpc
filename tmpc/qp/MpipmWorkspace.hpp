#pragma once

#include "detail/BundlePropertyMap.hpp"

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/core/Graph.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/Traits.hpp>

#include <tmpc/BlazeKernel.hpp>

#include <vector>
// #include <iostream>


namespace tmpc
{
    template <typename Real>
    class MpipmWorkspace
    {
    public:
        template <typename SizeMap>
        MpipmWorkspace(OcpGraph const& g, SizeMap size_map)
        :   graph_(g)
        ,   size_(num_vertices(g))
        {
            copyProperty(size_map, make_iterator_property_map(size_.begin(), vertexIndex(graph_)), vertices(graph_));

            // Allocate vertex properties of appropriate size
            vertexProperties_.reserve(num_vertices(graph_));
            for (auto const& sz : size_)
                vertexProperties_.emplace_back(sz);

            // Populate edgeIndex_ and allocate edge properties of appropriate size
            edgeProperties_.reserve(num_edges(graph_));
            for (auto e : edges(graph_))
                edgeProperties_.emplace_back(get(size_map, source(e, g)), get(size_map, target(e, g)));
        }


        auto size() const
        {
            return make_iterator_property_map(size_.begin(), vertexIndex(graph_));
        }


        auto const& graph() const
        {
            return graph_;
        }


        auto Q()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::Q_, vertexProperties());
        }


        auto R()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::R_, vertexProperties());
        }


        auto S()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::S_, vertexProperties());
        }


        auto q()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::q_, vertexProperties());
        }


        auto r()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::r_, vertexProperties());
        }


        auto lx()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lx_, vertexProperties());
        }


        auto ux()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ux_, vertexProperties());
        }


        auto lu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lu_, vertexProperties());
        }


        auto uu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::uu_, vertexProperties());
        }


        auto C()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::C_, vertexProperties());
        }


        auto D()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::D_, vertexProperties());
        }


        auto ld()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ld_, vertexProperties());
        }


        auto ud()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ud_, vertexProperties());
        }


        auto A()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::A_, edgeProperties());
        }


        auto B()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::B_, edgeProperties());
        }


        auto b()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::b_, edgeProperties());
        }


        auto x()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::x_, vertexProperties());
        }


        auto u()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::u_, vertexProperties());
        }


        auto lam_lx()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_lx_, vertexProperties());
        }


        auto lam_ux()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_ux_, vertexProperties());
        }


        auto lam_lu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_lu_, vertexProperties());
        }


        auto lam_uu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_uu_, vertexProperties());
        }


        auto lam_ld()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_ld_, vertexProperties());
        }


        auto lam_ud()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_ud_, vertexProperties());
        }


        auto pi()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::pi_, edgeProperties());
        }


        void solveUnconstrained()
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
                RiccatiBackwardVisitor(*this), 
                make_iterator_property_map(color.begin(), get(vertex_index, graph_)), 
                vertex(0, graph_));

            depth_first_search(graph_,
                RiccatiForwardVisitor(*this), 
                make_iterator_property_map(color.begin(), get(vertex_index, graph_)), 
                vertex(0, graph_));
        }


    private:
        struct VertexPropertyBundle
        {
            VertexPropertyBundle(OcpSize const& sz)
            :   Q_(sz.nx(), sz.nx())
            ,   R_(sz.nu(), sz.nu())
            ,   S_(sz.nu(), sz.nx())
            ,   q_(sz.nx())
            ,   r_(sz.nu())
            ,   lx_(sz.nx())
            ,   ux_(sz.nx())
            ,   lu_(sz.nu())
            ,   uu_(sz.nu())
            ,   C_(sz.nc(), sz.nx())
            ,   D_(sz.nc(), sz.nu())
            ,   ld_(sz.nc())
            ,   ud_(sz.nc())
            ,   x_(sz.nx())
            ,   u_(sz.nu())
            ,   lam_lx_(sz.nx())
            ,   lam_ux_(sz.nx())
            ,   lam_lu_(sz.nu())
            ,   lam_uu_(sz.nu())
            ,   lam_ld_(sz.nc())
            ,   lam_ud_(sz.nc())
            ,   P_(sz.nx(), sz.nx())
            ,   Lambda_(sz.nu(), sz.nu())
            ,   L_(sz.nu(), sz.nx())
            ,   p_(sz.nx())
            ,   l_(sz.nu())
            {
            }


            blaze::DynamicMatrix<Real> Q_;
            blaze::DynamicMatrix<Real> R_;
            blaze::DynamicMatrix<Real> S_;
            blaze::DynamicVector<Real> q_;
            blaze::DynamicVector<Real> r_;

            blaze::DynamicVector<Real> lx_;
            blaze::DynamicVector<Real> ux_;
            blaze::DynamicVector<Real> lu_;
            blaze::DynamicVector<Real> uu_;

            blaze::DynamicMatrix<Real> C_;
            blaze::DynamicMatrix<Real> D_;
            blaze::DynamicVector<Real> ld_;
            blaze::DynamicVector<Real> ud_;

            blaze::DynamicVector<Real> x_;
            blaze::DynamicVector<Real> u_;
            blaze::DynamicVector<Real> lam_lx_;
            blaze::DynamicVector<Real> lam_ux_;
            blaze::DynamicVector<Real> lam_lu_;
            blaze::DynamicVector<Real> lam_uu_;
            blaze::DynamicVector<Real> lam_ld_;
            blaze::DynamicVector<Real> lam_ud_;

            // For the Riccati solver
            blaze::DynamicMatrix<Real> P_;
            blaze::DynamicMatrix<Real> Lambda_;
            blaze::DynamicMatrix<Real> L_;
            blaze::DynamicVector<Real> p_;
            blaze::DynamicVector<Real> l_;
        };


        struct EdgePropertyBundle
        {
            EdgePropertyBundle(OcpSize const& size_src, OcpSize const& size_dst)
            :   A_(size_dst.nx(), size_src.nx())
            ,   B_(size_dst.nx(), size_src.nu())
            ,   b_(size_dst.nx())
            ,   pi_(size_dst.nx())
            {
            }


            blaze::DynamicMatrix<Real> A_;
            blaze::DynamicMatrix<Real> B_;
            blaze::DynamicVector<Real> b_;
            blaze::DynamicVector<Real> pi_;
        };


        class RiccatiBackwardVisitor
        :   public default_dfs_visitor 
        {
        public:
            RiccatiBackwardVisitor(MpipmWorkspace& ws)
            :   ws_(ws)
            {
            }
        
            
            void finish_vertex(OcpVertexDescriptor u, OcpGraph const& g) const
            {
                // std::clog << "Backward, Vertex " << u << std::endl;

                if (out_degree(u, g) == 0)
                {
                    // Alg 1 line 1
                    put(ws_.P(), u, get(ws_.Q(), u));

                    // Alg 2 line 1
                    put(ws_.p(), u, get(ws_.q(), u));
                }
                else
                {
                    auto const out_e = out_edges(u, g);

                    if (out_e.size() == 1)
                    {
                        auto const e = out_e.front();
                        auto const v = target(e, g);

                        // Alg 1 line 3
                        blaze::DynamicMatrix<Real> PA = trans(get(ws_.P(), v)) * get(ws_.A(), e);
                        blaze::DynamicMatrix<Real> PB = trans(get(ws_.P(), v)) * get(ws_.B(), e);

                        // Alg 1 line 4
                        blaze::DynamicMatrix<Real> BPA = trans(get(ws_.B(), e)) * PA;
                        blaze::DynamicMatrix<Real> BPB = trans(get(ws_.B(), e)) * PB;

                        // Alg 1 line 5
                        blaze::DynamicMatrix<Real> APA = trans(get(ws_.A(), e)) * PA;
                        
                        // Alg 1 line 6
                        // llh() or potrf()?
                        // TODO: llh() can be used with adaptors. See if using blaze::SymmetricMatrix improves the performance.
                        // https://bitbucket.org/blaze-lib/blaze/wiki/Matrix%20Operations#!cholesky-decomposition
                        //llh(get(ws_.Lambda(), u), get(ws_.Lambda(), u));
                        put(ws_.Lambda(), u, get(ws_.R(), u) + BPB);
                        potrf(get(ws_.Lambda(), u), 'L');

                        // Alg 1 line 7
                        put(ws_.L(), u, get(ws_.S(), u) + BPA);

                        if (!isEmpty(get(ws_.S(), u)))  // prevent trsm() from complaining about lda=0
                            trsm(get(ws_.Lambda(), u), get(ws_.L(), u), CblasLeft, CblasLower, 1.);
                        
                        // Alg 1 line 8
                        put(ws_.P(), u, get(ws_.Q(), u) + APA - trans(get(ws_.L(), u)) * get(ws_.L(), u));

                        // Alg 1 line 9
                        put(ws_.P(), u, 0.5 * (get(ws_.P(), u) + trans(get(ws_.P(), u))));

                        // Alg 2 line 3
                        put(ws_.l(), u, get(ws_.r(), u) + trans(get(ws_.B(), e)) * (trans(get(ws_.P(), v)) * get(ws_.b(), e) + get(ws_.p(), v)));
                        trsv(get(ws_.Lambda(), u), get(ws_.l(), u), 'L', 'N', 'N');

                        // Alg 2 line 4
                        put(ws_.p(), u, get(ws_.q(), u) + trans(get(ws_.A(), e)) * (trans(get(ws_.P(), v)) * get(ws_.b(), e) + get(ws_.p(), v)) 
                            - trans(get(ws_.L(), u)) * get(ws_.l(), u));
                    }
                    else
                    {
                        throw std::invalid_argument("Riccati solver is not implemented on tree QPs yet");
                    }
                }

                // std::clog << "Lambda = " << std::endl << get(ws_.Lambda(), u) << std::endl;
                // std::clog << "L = " << std::endl << get(ws_.L(), u) << std::endl;
                // std::clog << "l = " << std::endl << get(ws_.l(), u) << std::endl;
                // std::clog << "P = " << std::endl << get(ws_.P(), u) << std::endl;
                // std::clog << "p = " << std::endl << get(ws_.p(), u) << std::endl;
            }


        private:
            MpipmWorkspace& ws_;
        };


        class RiccatiForwardVisitor
        :   public default_dfs_visitor 
        {
        public:
            RiccatiForwardVisitor(MpipmWorkspace& ws)
            :   ws_(ws)
            {
            }
        
            
            void discover_vertex(OcpVertexDescriptor u, OcpGraph const& g) const
            {
                // std::clog << "Forward, Vertex " << u << std::endl;

                if (in_degree(u, g) == 0)
                {
                    // Root vertex
                    put(ws_.x(), u, -inv(get(ws_.P(), u)) * get(ws_.p(), u));
                }

                // Alg 2 line 8
                put(ws_.u(), u, -(get(ws_.L(), u) * get(ws_.x(), u) + get(ws_.l(), u)));
                trsv(get(ws_.Lambda(), u), get(ws_.u(), u), 'L', 'T', 'N');

                /*
                std::clog << "u = " << std::endl << get(ws_.u(), u) << std::endl;
                std::clog << "x = " << std::endl << get(ws_.x(), u) << std::endl;
                */
            }


            void tree_edge(OcpEdgeDescriptor e, OcpGraph const& g) const
            {
                // std::clog << "Forward, Edge " << e << std::endl;
                auto const u = source(e, g);
                auto const v = target(e, g);

                // Alg 2 line 9
                put(ws_.x(), v, get(ws_.A(), e) * get(ws_.x(), u) + get(ws_.B(), e) * get(ws_.u(), u) + get(ws_.b(), e));

                // Alg 2 line 10
                put(ws_.pi(), e, trans(get(ws_.P(), v)) * get(ws_.x(), v) + get(ws_.p(), v));

                // std::clog << "pi = " << std::endl << get(ws_.pi(), e) << std::endl;
            }


        private:
            MpipmWorkspace& ws_;
        };


        auto vertexProperties()
        {
            return make_iterator_property_map(vertexProperties_.begin(), vertexIndex(graph_));
        }


        auto edgeProperties()
        {
            return make_iterator_property_map(edgeProperties_.begin(), edgeIndex());
        }


        auto edgeIndex() const
        {
            return get(edge_index, graph_);
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
        std::vector<EdgePropertyBundle> edgeProperties_;
    };


    template <typename Real>
    struct KernelOf<MpipmWorkspace<Real>>
    {
        using type = BlazeKernel<Real>;
    };


    template <typename Real>
    struct RealOf<MpipmWorkspace<Real>>
    {
        using type = Real;
    };
}