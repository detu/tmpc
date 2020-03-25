#pragma once

#include <tmpc/property_map/BundlePropertyMap.hpp>

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
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


        auto H()
        {
            return BundlePropertyMap(&VertexPropertyBundle::H_, vertexProperties());
        }


        auto H() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::H_, vertexProperties());
        }


        auto Q()
        {
            return BundlePropertyMap(&VertexPropertyBundle::Q_, vertexProperties());
        }


        auto Q() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::Q_, vertexProperties());
        }


        auto R()
        {
            return BundlePropertyMap(&VertexPropertyBundle::R_, vertexProperties());
        }


        auto R() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::R_, vertexProperties());
        }


        auto S()
        {
            return BundlePropertyMap(&VertexPropertyBundle::S_, vertexProperties());
        }


        auto S() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::S_, vertexProperties());
        }


        auto q()
        {
            return BundlePropertyMap(&VertexPropertyBundle::q_, vertexProperties());
        }


        auto q() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::q_, vertexProperties());
        }


        auto r()
        {
            return BundlePropertyMap(&VertexPropertyBundle::r_, vertexProperties());
        }


        auto r() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::r_, vertexProperties());
        }


        auto lx()
        {
            return BundlePropertyMap(&VertexPropertyBundle::lx_, vertexProperties());
        }


        auto lx() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::lx_, vertexProperties());
        }


        auto ux()
        {
            return BundlePropertyMap(&VertexPropertyBundle::ux_, vertexProperties());
        }


        auto ux() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::ux_, vertexProperties());
        }


        auto lu()
        {
            return BundlePropertyMap(&VertexPropertyBundle::lu_, vertexProperties());
        }


        auto lu() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::lu_, vertexProperties());
        }


        auto uu()
        {
            return BundlePropertyMap(&VertexPropertyBundle::uu_, vertexProperties());
        }


        auto uu() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::uu_, vertexProperties());
        }


        auto C()
        {
            return BundlePropertyMap(&VertexPropertyBundle::C_, vertexProperties());
        }


        auto C() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::C_, vertexProperties());
        }


        auto D()
        {
            return BundlePropertyMap(&VertexPropertyBundle::D_, vertexProperties());
        }


        auto D() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::D_, vertexProperties());
        }


        auto ld()
        {
            return BundlePropertyMap(&VertexPropertyBundle::ld_, vertexProperties());
        }


        auto ld() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::ld_, vertexProperties());
        }


        auto ud()
        {
            return BundlePropertyMap(&VertexPropertyBundle::ud_, vertexProperties());
        }


        auto ud() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::ud_, vertexProperties());
        }



		/*
        auto BA()
        {
            return BundlePropertyMap(&EdgePropertyBundle::BA_, edgeProperties());
        }


        auto BA() const
        {
            return BundlePropertyMap(&EdgePropertyBundle::BA_, edgeProperties());
        }
		*/


        auto A()
        {
            return BundlePropertyMap(&EdgePropertyBundle::A_, edgeProperties());
        }


        auto A() const
        {
            return BundlePropertyMap(&EdgePropertyBundle::A_, edgeProperties());
        }


        auto B()
        {
            return BundlePropertyMap(&EdgePropertyBundle::B_, edgeProperties());
        }


        auto B() const
        {
            return BundlePropertyMap(&EdgePropertyBundle::B_, edgeProperties());
        }


        auto b()
        {
            return BundlePropertyMap(&EdgePropertyBundle::b_, edgeProperties());
        }


        auto b() const
        {
            return BundlePropertyMap(&EdgePropertyBundle::b_, edgeProperties());
        }


        auto x()
        {
            return BundlePropertyMap(&VertexPropertyBundle::x_, vertexProperties());
        }


        auto x() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::x_, vertexProperties());
        }


        auto u()
        {
            return BundlePropertyMap(&VertexPropertyBundle::u_, vertexProperties());
        }


        auto u() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::u_, vertexProperties());
        }


        auto lam_lx()
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_lx_, vertexProperties());
        }


        auto lam_lx() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_lx_, vertexProperties());
        }


        auto lam_ux()
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_ux_, vertexProperties());
        }


        auto lam_ux() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_ux_, vertexProperties());
        }


        auto lam_lu()
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_lu_, vertexProperties());
        }


        auto lam_lu() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_lu_, vertexProperties());
        }


        auto lam_uu()
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_uu_, vertexProperties());
        }


        auto lam_uu() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_uu_, vertexProperties());
        }


        auto lam_ld()
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_ld_, vertexProperties());
        }


        auto lam_ld() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_ld_, vertexProperties());
        }


        auto lam_ud()
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_ud_, vertexProperties());
        }


        auto lam_ud() const
        {
            return BundlePropertyMap(&VertexPropertyBundle::lam_ud_, vertexProperties());
        }


        auto pi()
        {
            return BundlePropertyMap(&EdgePropertyBundle::pi_, edgeProperties());
        }


        auto pi() const
        {
            return BundlePropertyMap(&EdgePropertyBundle::pi_, edgeProperties());
        }


    private:
        struct VertexPropertyBundle
        {
            VertexPropertyBundle(OcpSize const& sz)
            :   H_(sz.nu() + sz.nx())
            ,   Q_(submatrix(H_, sz.nu(), sz.nu(), sz.nx(), sz.nx()))
            ,   R_(submatrix(H_, 0, 0, sz.nu(), sz.nu()))
            ,   S_(submatrix(H_, 0, sz.nu(), sz.nu(), sz.nx()))
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

            // H = [ R,   S
            //       S^T, Q]
            blaze::SymmetricMatrix<blaze::DynamicMatrix<Real, blaze::columnMajor>> H_;

            using Submatrix = decltype(blaze::submatrix(H_, 0, 0, 1, 1));

            Submatrix Q_;
            Submatrix R_;
            Submatrix S_;
            blaze::DynamicVector<Real> q_;
            blaze::DynamicVector<Real> r_;

            blaze::DynamicVector<Real> lx_;
            blaze::DynamicVector<Real> ux_;
            blaze::DynamicVector<Real> lu_;
            blaze::DynamicVector<Real> uu_;

            blaze::DynamicMatrix<Real, blaze::columnMajor> C_;
            blaze::DynamicMatrix<Real, blaze::columnMajor> D_;
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
			/*
            :   BA_(size_dst.nx(), size_src.nu() + size_src.nx())
            ,   A_(submatrix(BA_, 0, size_src.nu(), size_dst.nx(), size_src.nx()))
            ,   B_(submatrix(BA_, 0, 0, size_dst.nx(), size_src.nu()))
			*/
            ,   b_(size_dst.nx())
            ,   pi_(size_dst.nx())
            {
            }


            blaze::DynamicMatrix<Real, blaze::columnMajor> A_;
            blaze::DynamicMatrix<Real, blaze::columnMajor> B_;
			
			/*
            blaze::DynamicMatrix<Real, blaze::columnMajor> BA_;

            using Submatrix = decltype(blaze::submatrix(BA_, 0, 0, 1, 1));

            Submatrix A_;
            Submatrix B_;
			*/
			
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