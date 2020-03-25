#pragma once

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/property_map/BundlePropertyMap.hpp>

#include <tmpc/Math.hpp>

#include <vector>


namespace tmpc
{
    template <typename Real>
    class OcpKktResidual
    {
    public:
        template <typename SizeMap>
        OcpKktResidual(OcpGraph const& g, SizeMap size_map)
        :   graph_ {g}
        ,   vertexData_(num_vertices(g))
        ,   edgeData_(num_edges(g))
        {
            for (auto v : graph::vertices(g))
                vertexData_[get(graph::vertex_index, g, v)].resize(get(size_map, v));

            for (auto e : graph::edges(g))
                edgeData_[get(graph::edge_index, g, e)].resize(
                    get(size_map, source(e, g)), get(size_map, target(e, g)));
        }


        auto gx() const noexcept
        {
            return BundlePropertyMap(&VertexData::gx_, make_iterator_property_map(begin(vertexData_), get(graph::vertex_index, graph_)));
        }


        auto gx()
        {
            return BundlePropertyMap(&VertexData::gx_, make_iterator_property_map(begin(vertexData_), get(graph::vertex_index, graph_)));
        }


        auto gu() const noexcept
        {
            return BundlePropertyMap(&VertexData::gu_, make_iterator_property_map(begin(vertexData_), get(graph::vertex_index, graph_)));
        }


        auto gu()
        {
            return BundlePropertyMap(&VertexData::gu_, make_iterator_property_map(begin(vertexData_), get(graph::vertex_index, graph_)));
        }


        auto c() const noexcept
        {
            return BundlePropertyMap(&EdgeData::c_, make_iterator_property_map(begin(edgeData_), get(graph::edge_index, graph_)));
        }


        auto c()
        {
            return BundlePropertyMap(&EdgeData::c_, make_iterator_property_map(begin(edgeData_), get(graph::edge_index, graph_)));
        }
        

    private:
        struct VertexData
        {
            VertexData() = default;


            VertexData(OcpSize const& sz)
            :   gx_ {sz.nx()}
            ,   gu_ {sz.nu()}
            {
            }


            void resize(OcpSize const& sz)
            {
                gx_.resize(sz.nx());
                gu_.resize(sz.nu());
            }


            // Gradient of the Lagrangian w.r.t. x
            blaze::DynamicVector<Real, blaze::columnVector> gx_;

            // Gradient of the Lagrangian w.r.t. u
            blaze::DynamicVector<Real, blaze::columnVector> gu_;
        };


        struct EdgeData
        {
            EdgeData() = default;


            EdgeData(OcpSize const& sz_src, OcpSize const& sz_dst)
            :   c_ {sz_dst.nx()}
            {
            }


            void resize(OcpSize const& sz_src, OcpSize const& sz_dst)
            {
                c_.resize(sz_dst.nx());
            }


            // Equality constraint residual
            blaze::DynamicVector<Real, blaze::columnVector> c_;
        };


        OcpGraph graph_;

        std::vector<VertexData> vertexData_;
        std::vector<EdgeData> edgeData_;
    };
}