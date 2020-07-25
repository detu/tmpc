/// @brief Demonstrates working with OCP trees:
/// creating, traversal, vertex and edge functions, copying, comparison.

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/print/ocp/OcpTree.hpp>

#include <iostream>
#include <array>


int main(int, char **)
{
    using namespace tmpc;

    OcpTree g {2,2,2,1,1,1,1,1,1,1,1,0,0,0,0};

    for (auto v : vertices(g))
    {
        std::cout << "Vertex " << v;

        auto const kids = children(v, g);
        if (std::empty(kids))
            std::cout << " has no children";
        else
            std::cout << " is the parent of";

        for (auto k : kids) 
            std::cout << " " << k;

        std::cout << " and is at depth " << g.depth(v) << std::endl;
    }

    for (auto e : edges(g))
        std::cout << "Edge " << e << " from " << source(e, g) << " to " << target(e, g) << std::endl;

    OcpTree const g1 = g, g2 {2, 1, 1, 0, 0};
    std::cout << "g1 == g: " << (g1 == g ? "yes" : "no") << std::endl;
    std::cout << "g2 == g: " << (g2 == g ? "yes" : "no") << std::endl;
    std::cout << "g1:" << std::endl << g1;

    return EXIT_SUCCESS;
}