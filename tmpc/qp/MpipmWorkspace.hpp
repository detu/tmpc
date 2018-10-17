#pragma once

#include "detail/BundlePropertyMap.hpp"

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/PropertyMap.hpp>
#include <tmpc/graph/Graph.hpp>
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
            copyProperty(size_map, make_iterator_property_map(size_.begin(), vertexIndex(graph_)), graph::vertices(graph_));

            // Allocate vertex properties of appropriate size
            vertexProperties_.reserve(num_vertices(graph_));
            for (auto const& sz : size_)
                vertexProperties_.emplace_back(sz);

            // Populate edgeIndex_ and allocate edge properties of appropriate size
            edgeProperties_.reserve(num_edges(graph_));
            for (auto e : graph::edges(graph_))
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


        auto Q() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::Q_, vertexProperties());
        }


        auto R()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::R_, vertexProperties());
        }


        auto R() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::R_, vertexProperties());
        }


        auto S()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::S_, vertexProperties());
        }


        auto S() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::S_, vertexProperties());
        }


        auto q()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::q_, vertexProperties());
        }


        auto q() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::q_, vertexProperties());
        }


        auto r()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::r_, vertexProperties());
        }


        auto r() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::r_, vertexProperties());
        }


        auto lx()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lx_, vertexProperties());
        }


        auto lx() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lx_, vertexProperties());
        }


        auto ux()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ux_, vertexProperties());
        }


        auto ux() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ux_, vertexProperties());
        }


        auto lu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lu_, vertexProperties());
        }


        auto lu() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lu_, vertexProperties());
        }


        auto uu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::uu_, vertexProperties());
        }


        auto uu() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::uu_, vertexProperties());
        }


        auto C()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::C_, vertexProperties());
        }


        auto C() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::C_, vertexProperties());
        }


        auto D()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::D_, vertexProperties());
        }


        auto D() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::D_, vertexProperties());
        }


        auto ld()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ld_, vertexProperties());
        }


        auto ld() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ld_, vertexProperties());
        }


        auto ud()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ud_, vertexProperties());
        }


        auto ud() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::ud_, vertexProperties());
        }


        auto A()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::A_, edgeProperties());
        }


        auto A() const
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::A_, edgeProperties());
        }


        auto B()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::B_, edgeProperties());
        }


        auto B() const
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::B_, edgeProperties());
        }


        auto b()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::b_, edgeProperties());
        }


        auto b() const
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::b_, edgeProperties());
        }


        auto x()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::x_, vertexProperties());
        }


        auto x() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::x_, vertexProperties());
        }


        auto u()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::u_, vertexProperties());
        }


        auto u() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::u_, vertexProperties());
        }


        auto lam_lx()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_lx_, vertexProperties());
        }


        auto lam_lx() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_lx_, vertexProperties());
        }


        auto lam_ux()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_ux_, vertexProperties());
        }


        auto lam_ux() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_ux_, vertexProperties());
        }


        auto lam_lu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_lu_, vertexProperties());
        }


        auto lam_lu() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_lu_, vertexProperties());
        }


        auto lam_uu()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_uu_, vertexProperties());
        }


        auto lam_uu() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_uu_, vertexProperties());
        }


        auto lam_ld()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_ld_, vertexProperties());
        }


        auto lam_ld() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_ld_, vertexProperties());
        }


        auto lam_ud()
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_ud_, vertexProperties());
        }


        auto lam_ud() const
        {
            return detail::BundlePropertyMap(&VertexPropertyBundle::lam_ud_, vertexProperties());
        }


        auto pi()
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::pi_, edgeProperties());
        }


        auto pi() const
        {
            return detail::BundlePropertyMap(&EdgePropertyBundle::pi_, edgeProperties());
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


        auto vertexProperties()
        {
            return make_iterator_property_map(vertexProperties_.begin(), vertexIndex(graph_));
        }


        auto vertexProperties() const
        {
            return make_iterator_property_map(vertexProperties_.begin(), vertexIndex(graph_));
        }


        auto edgeProperties()
        {
            return make_iterator_property_map(edgeProperties_.begin(), edgeIndex());
        }


        auto edgeProperties() const
        {
            return make_iterator_property_map(edgeProperties_.begin(), edgeIndex());
        }


        auto edgeIndex() const
        {
            return get(graph::edge_index, graph_);
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