#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/wavefront.hpp>

#include <tmpc/ocp/OcpSizeGraph.hpp>
#include <tmpc/core/IteratorRange.hpp>
#include <tmpc/core/GraphTools.hpp>


int main()
{
    using namespace boost;
    using namespace tmpc;

    size_t const N = 4;
    OcpSizeGraph g(N);
    add_edge(0, 1, 0, g);
    add_edge(1, 2, 1, g);
    add_edge(0, 3, 2, g);

    auto index_map = get(vertex_index, g);

    for (auto i : verticesR(g))
    {
        std::cout << get(index_map, i);

        auto adj_vert = adjacent_verticesR(i, g);
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

    g[0].size = OcpSize {2, 3, 1, 1};
    g[2].size = OcpSize {3, 4, 2, 3};

    for (auto v : verticesR(g))
    {
        auto const& sz = g[v].size;
        std::cout << "Vertex " << v << ", index=" << get(index_map, v)
            << ", nx=" << sz.nx() << ", nu=" << sz.nu() << ", nc=" << sz.nc() << ", ns=" << sz.ns() << std::endl;
    }

    std::cout << "Max wavefront = " << max_wavefront(g) << std::endl;

    auto edge_index_map = get(edge_index, g);
    for (auto e : edgesR(g))
        std::cout << "Edge " << e << ", index=" << get(edge_index_map, e) << std::endl;

    remove_vertex(0, g);
    std::cout << "After removing vertex 0:" << std::endl;
    for (auto v : verticesR(g))
    {
        auto const& sz = g[v].size;
        std::cout << "Vertex " << v << ", index=" << get(index_map, v)
            << ", nx=" << sz.nx() << ", nu=" << sz.nu() << ", nc=" << sz.nc() << ", ns=" << sz.ns() << std::endl;
    }

    for (auto e : edgesR(g))
        std::cout << "Edge " << e << ", index=" << get(edge_index_map, e) << std::endl;

    return EXIT_SUCCESS;
}