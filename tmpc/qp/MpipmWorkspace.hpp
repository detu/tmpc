#pragma once

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/core/PropertyMap.hpp>

#include <blaze/Math.h>

#include <vector>


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
            copyProperty(size_map, make_iterator_property_map(size_.begin(), vertexIndex(g)), vertices(g));
        }


    private:
        OcpGraph graph_;
        std::vector<OcpSize> size_;
    };
}