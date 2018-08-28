#include <tmpc/ocp/OcpSizeGraph.hpp>


namespace tmpc 
{
	OcpSizeGraph ocpSizeGraphNominalMpc(size_t nt, size_t nx, size_t nu, size_t nc, size_t nct)
	{
		// Create the graph.
        OcpSizeGraph g(nt + 1);

		// All but the last stages.
		size_t v = 0;
        for (; v < nt; ++v)
        {
            g[v].size = OcpSize {nx, nu, nc};

            if (v > 0)
                add_edge(v - 1, v, v - 1 /* edge index */, g);
        }

		// Last stage.
		g[v].size = OcpSize {nx, nu, nct};
		if (v > 0)
            add_edge(v - 1, v, v - 1 /* edge index */, g);

        return g;
	}
}
