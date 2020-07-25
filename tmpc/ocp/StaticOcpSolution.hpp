#pragma once

#include <tmpc/property_map/BundlePropertyMap.hpp>

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/StaticOcpSize.hpp>
#include <tmpc/property_map/PropertyMap.hpp>

#include <tmpc/Math.hpp>

#include <vector>


namespace tmpc
{
    /// @brief Solution of an OCP with sizes known at compile-time.
    ///
    template <typename R, size_t NX, size_t NU, size_t NC = 0>
    class StaticOcpSolution
    {
    public:
        using Real = R;

        
        StaticOcpSolution(OcpTree const& g)
        :   graph_ {g}
        ,   size_ {g}
        ,   vertexProperties_(num_vertices(g))
        ,   edgeProperties_(num_edges(g))
        {
        }


        auto const& size() const
        {
            return size_;
        }


        auto const& graph() const
        {
            return graph_;
        }


        auto const& x(OcpVertex v) const
        {
            return vertexProperties_[v].x_;
        }


        auto& x(OcpVertex v)
        {
            return vertexProperties_[v].x_;
        }


        auto const& u(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].u_;
        }


        auto& u(OcpVertex v)
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].u_;
        }


        auto const& lam_lx(OcpVertex v) const
        {
            return vertexProperties_[v].lam_lx_;
        }


        auto& lam_lx(OcpVertex v)
        {
            return vertexProperties_[v].lam_lx_;
        }


        auto const& lam_ux(OcpVertex v) const
        {
            return vertexProperties_[v].lam_ux_;
        }


        auto& lam_ux(OcpVertex v)
        {
            return vertexProperties_[v].lam_ux_;
        }


        auto const& lam_lu(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].lam_lu_;
        }


        auto& lam_lu(OcpVertex v)
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].lam_lu_;
        }


        auto const& lam_uu(OcpVertex v) const
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].lam_uu_;
        }


        auto& lam_uu(OcpVertex v)
        {
            assert(out_degree(v, graph_) > 0);
            return vertexProperties_[v].lam_uu_;
        }


        auto const& lam_ld(OcpVertex v) const
        {
            return vertexProperties_[v].lam_ld_;
        }


        auto& lam_ld(OcpVertex v)
        {
            return vertexProperties_[v].lam_ld_;
        }


        auto const& lam_ud(OcpVertex v) const
        {
            return vertexProperties_[v].lam_ud_;
        }


        auto& lam_ud(OcpVertex v)
        {
            return vertexProperties_[v].lam_ud_;
        }


        auto const& pi(OcpEdge e) const
        {
            return edgeProperties_[e].pi_;
        }


        auto& pi(OcpEdge e)
        {
            return edgeProperties_[e].pi_;
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


        OcpTree graph_;
        StaticOcpSize<NX, NU, NC> size_;
        std::vector<VertexPropertyBundle> vertexProperties_;
        std::vector<EdgePropertyBundle> edgeProperties_;
    };
}