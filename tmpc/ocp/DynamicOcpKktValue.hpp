#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/property_map/BundlePropertyMap.hpp>

#include <tmpc/Math.hpp>

#include <vector>


namespace tmpc
{
    template <typename R>
    class DynamicOcpKktValue
    {
    public:
        using Real = R;
        

        template <OcpSize Size>
        DynamicOcpKktValue(Size const& size)
        :   size_ {size}
        ,   vertexData_(num_vertices(graph()))
        ,   edgeData_(num_edges(graph()))
        {
            for (auto v : vertices(graph()))
                vertexData_[v].resize(size, v);

            for (auto e : edges(graph()))
                edgeData_[e].resize(size, e);
        }


        auto const& graph() const noexcept
        {
            return size_.graph();
        }


        auto const& size() const noexcept
        {
            return size_;
        }


        auto const& gx(OcpVertex v) const
        {
            return vertexData_[v].gx_;
        }


        template <typename VT>
        void gx(OcpVertex v, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            noresize(vertexData_[v].gx_) = ~val;
        }


        auto const& gu(OcpVertex v) const
        {
            return vertexData_[v].gu_;
        }


        template <typename VT>
        void gu(OcpVertex v, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            noresize(vertexData_[v].gu_) = ~val;
        }


        auto const& c(OcpEdge e) const
        {
            return edgeData_[e].c_;
        }


        template <typename VT>
        void c(OcpEdge e, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            noresize(edgeData_[e].c_) = ~val;
        }
        

    private:
        struct VertexData
        {
            VertexData() = default;


            VertexData(DynamicOcpSize const& sz, OcpVertex v)
            :   gx_ {sz.nx(v)}
            ,   gu_ {sz.nu(v)}
            {
            }


            void resize(DynamicOcpSize const& sz, OcpVertex v)
            {
                gx_.resize(sz.nx(v));
                gu_.resize(sz.nu(v));
            }


            // Gradient of the Lagrangian w.r.t. x
            blaze::DynamicVector<Real, blaze::columnVector> gx_;

            // Gradient of the Lagrangian w.r.t. u
            blaze::DynamicVector<Real, blaze::columnVector> gu_;
        };


        struct EdgeData
        {
            EdgeData() = default;


            EdgeData(DynamicOcpSize const& sz, OcpEdge e)
            :   c_ {sz.nx(target(e, sz.graph()))}
            {
            }


            void resize(DynamicOcpSize const& sz, OcpEdge e)
            {
                c_.resize(sz.nx(target(e, sz.graph())));
            }


            // Equality constraint residual
            blaze::DynamicVector<Real, blaze::columnVector> c_;
        };


        DynamicOcpSize size_;

        std::vector<VertexData> vertexData_;
        std::vector<EdgeData> edgeData_;
    };
}