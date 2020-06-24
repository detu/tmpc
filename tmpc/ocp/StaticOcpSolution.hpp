#pragma once

#include <tmpc/property_map/BundlePropertyMap.hpp>

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/graph/Graph.hpp>
#include <tmpc/Math.hpp>

#include <vector>


namespace tmpc
{
    /// @brief Solution of an OCP with sizes known at compile-time.
    ///
    template <typename Real, size_t NX, size_t NU, size_t NC = 0>
    class StaticOcpSolution
    {
    public:
        StaticOcpSolution(OcpGraph const& g)
        :   graph_(g)
        ,   vertexProperties_(num_vertices(g))
        ,   edgeProperties_(num_edges(g))
        {
        }


        auto size() const
        {
            return make_function_property_map<OcpVertexDescriptor>(
                [] (OcpVertexDescriptor v) 
                { 
                    return OcpSize {NX, NU, NC};
                }
            );
        }


        auto const& graph() const
        {
            return graph_;
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
            blaze::StaticVector<Real, NX> x_;
            blaze::StaticVector<Real, NU> u_;
            blaze::StaticVector<Real, NX> lam_lx_;
            blaze::StaticVector<Real, NX> lam_ux_;
            blaze::StaticVector<Real, NU> lam_lu_;
            blaze::StaticVector<Real, NU> lam_uu_;
            blaze::StaticVector<Real, NC> lam_ld_;
            blaze::StaticVector<Real, NC> lam_ud_;
        };


        struct EdgePropertyBundle
        {
            blaze::StaticVector<Real, NX> b_;
            blaze::StaticVector<Real, NX> pi_;
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
        std::vector<VertexPropertyBundle> vertexProperties_;
        std::vector<EdgePropertyBundle> edgeProperties_;
    };
}