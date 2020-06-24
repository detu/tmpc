#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/wavefront.hpp>

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/property_map/PropertyMap.hpp>


int main()
{
    using namespace tmpc;

    size_t const N = 4;
    std::vector<std::pair<size_t, size_t>> edge_list {{0, 1}, {1, 2}, {0, 3}};
    OcpGraph g(edge_list.begin(), edge_list.end(), N);

    auto index_map = vertexIndex(g);

    for (auto i : graph::vertices(g))
    {
        std::cout << get(index_map, i);

        auto adj_vert = graph::adjacent_vertices(i, g);
        if (adj_vert.size() == 0)
            std::cout << " has no children";
        else
            std::cout << " is the parent of ";

        for (auto ai = adj_vert.begin(); ai != adj_vert.end(); ++ai) 
        {
            std::cout << get(index_map, *ai);
            if (std::next(ai) != adj_vert.end())
                std::cout << ", ";
        }

        std::cout << std::endl;
    }

    std::vector<OcpSize> const size = 
    {
        OcpSize {2, 3, 1, 1},
        OcpSize {3, 4, 2, 3}
    };

    auto const size_map = iterator_property_map(size.begin(), index_map);

    for (auto v : graph::vertices(g))
    {
        auto const& sz = size_map[v];
        std::cout << "Vertex " << v << ", index=" << get(index_map, v)
            << ", nx=" << sz.nx() << ", nu=" << sz.nu() << ", nc=" << sz.nc() << ", ns=" << sz.ns() << std::endl;
    }

    std::cout << "Max wavefront = " << max_wavefront(g) << std::endl;

    for (auto e : graph::edges(g))
        std::cout << "Edge " << get(graph::edge_index, g, e) << std::endl;

    return EXIT_SUCCESS;
}