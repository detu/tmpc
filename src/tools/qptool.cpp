#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>

#include <nlohmann/json.hpp>

#include <tmpc/ocp/OcpSizeGraph.hpp>
#include <tmpc/core/IteratorRange.hpp>

#include <boost/range/adaptor/indexed.hpp>


using namespace tmpc;


void printSizeGraph(OcpSizeGraph const& g)
{
    auto index_map = get(boost::vertex_index, g);
    auto ocp_size = get(OcpSizeProperty_t(), g);

    for (auto i : make_iterator_range(vertices(g)))
    {
        std::cout << get(index_map, i);
        auto adj_vert = make_iterator_range(adjacent_vertices(i, g));
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

    for (auto v : make_iterator_range(vertices(g)))
    {
        auto sz = get(ocp_size, v);
        std::cout << "Vertex " << v << ", nx=" << sz.nx() << ", nu=" << sz.nu() << ", nc=" << sz.nc() << ", ns=" << sz.ns() << std::endl;
    }
}


int main(int argc, char ** argv)
{
    using nlohmann::json;
    using boost::adaptors::indexed;

    json j;
    std::cin >> j;

    size_t const N = j["nodes"].size();
    OcpSizeGraph g(N);

    auto ocp_size = get(OcpSizeProperty_t(), g);
    for (auto const& j_vertex : j["nodes"] | indexed(0))
    {
        auto const& j_v = j_vertex.value();
        put(ocp_size, j_vertex.index(), OcpSize {j_v["q"].size(), j_v["r"].size(), j_v["ld"].size(), j_v["zl"].size()});
    }

    for (auto j_edge : j["edges"])
        add_edge(j_edge["from"], j_edge["to"], g);

    printSizeGraph(g);

    return EXIT_SUCCESS;
}