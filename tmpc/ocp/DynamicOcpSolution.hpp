#pragma once

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/property_map/BundlePropertyMap.hpp>

#include <vector>


namespace tmpc
{
    template <typename RT>
    class DynamicOcpSolution
    {
    public:
        using Real = RT;


        template <OcpSize Size>
        DynamicOcpSolution(Size const& size)
        :   graph_ {size.graph()}
        ,   size_ {size}
        {
            // Allocate vertex properties of appropriate size
            vertexProperties_.reserve(num_vertices(graph_));
            for (auto v : vertices(graph_))
                vertexProperties_.emplace_back(size_, v);

            // Populate edgeIndex_ and allocate edge properties of appropriate size
            edgeProperties_.reserve(num_edges(graph_));
            for (auto e : edges(graph_))
                edgeProperties_.emplace_back(size_, e);
        }


        auto const& size() const noexcept
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


        auto x(OcpVertex v)
        {
            return noresize(vertexProperties_[v].x_);
        }


        auto const& u(OcpVertex v) const
        {
            return vertexProperties_[v].u_;
        }


        auto u(OcpVertex v)
        {
            return noresize(vertexProperties_[v].u_);
        }


        auto const& lam_lx(OcpVertex v) const
        {
            return vertexProperties_[v].lam_lx_;
        }


        auto lam_lx(OcpVertex v)
        {
            return noresize(vertexProperties_[v].lam_lx_);
        }


        auto const& lam_ux(OcpVertex v) const
        {
            return vertexProperties_[v].lam_ux_;
        }


        auto lam_ux(OcpVertex v)
        {
            return noresize(vertexProperties_[v].lam_ux_);
        }


        auto const& lam_lu(OcpVertex v) const
        {
            return vertexProperties_[v].lam_lu_;
        }


        auto lam_lu(OcpVertex v)
        {
            return noresize(vertexProperties_[v].lam_lu_);
        }


        auto const& lam_uu(OcpVertex v) const
        {
            return vertexProperties_[v].lam_uu_;
        }


        auto lam_uu(OcpVertex v)
        {
            return noresize(vertexProperties_[v].lam_uu_);
        }


        auto const& lam_ld(OcpVertex v) const
        {
            return vertexProperties_[v].lam_ld_;
        }


        auto lam_ld(OcpVertex v)
        {
            return noresize(vertexProperties_[v].lam_ld_);
        }


        auto const& lam_ud(OcpVertex v) const
        {
            return vertexProperties_[v].lam_ud_;
        }


        auto lam_ud(OcpVertex v)
        {
            return noresize(vertexProperties_[v].lam_ud_);
        }


        auto const& pi(OcpEdge e) const
        {
            return edgeProperties_[e].pi_;
        }


        auto pi(OcpEdge e)
        {
            return noresize(edgeProperties_[e].pi_);
        }


    private:
        struct VertexPropertyBundle
        {
            VertexPropertyBundle(DynamicOcpSize const& sz, OcpVertex v)
            :   x_(sz.nx(v))
            ,   u_(sz.nu(v))
            ,   lam_lx_(sz.nx(v))
            ,   lam_ux_(sz.nx(v))
            ,   lam_lu_(sz.nu(v))
            ,   lam_uu_(sz.nu(v))
            ,   lam_ld_(sz.nc(v))
            ,   lam_ud_(sz.nc(v))
            {
            }


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
            EdgePropertyBundle(DynamicOcpSize const& size, OcpEdge e)
            :   pi_(size.nx(target(e, size.graph())))
            {
            }


            blaze::DynamicVector<Real> pi_;
        };


        OcpTree graph_;
        DynamicOcpSize size_;
        std::vector<VertexPropertyBundle> vertexProperties_;
        std::vector<EdgePropertyBundle> edgeProperties_;
    };
}