#include <tmpc/ocp/OcpGraph.hpp>


namespace tmpc
{
    OcpGraph ocpGraphLinear(size_t n)
    {
        // Create the graph.
        OcpGraph g(n);

        for (size_t v = 0; v < n; ++v)
        {
            if (v > 0)
                add_edge(v - 1, v, g);
        }

        return g;
    }
}