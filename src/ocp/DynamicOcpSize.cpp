#include <tmpc/ocp/DynamicOcpSize.hpp>


namespace tmpc
{
    DynamicOcpSize::DynamicOcpSize(std::initializer_list<OcpVertexSize> vertex_size)
    :	graph_(vertex_size.size())
    ,	size_ {new OcpVertexSize[vertex_size.size()]}
    {
        std::copy(begin(vertex_size), end(vertex_size), size_.get());
    }


    DynamicOcpSize::DynamicOcpSize(OcpTree const& graph, std::initializer_list<OcpVertexSize> vertex_size)
    :	graph_(graph)
    ,	size_ {new OcpVertexSize[vertex_size.size()]}
    {
        if (vertex_size.size() != num_vertices(graph_))
            TMPC_THROW_EXCEPTION(std::invalid_argument(
                "The size of vertex_size must be equal to number of vertices"));

        std::copy(begin(vertex_size), end(vertex_size), size_.get());
    }


    DynamicOcpSize::DynamicOcpSize(OcpTree const& graph, size_t nx, size_t nu, 
        size_t nc, size_t ns, size_t nct, bool first_state_empty)
    :	DynamicOcpSize {
            graph,
            [nx, nu, nc, ns, nct, first_state_empty] (OcpVertex v, OcpTree const& g)
            {
                size_t const actual_nx = v == root(g) && first_state_empty ? 0 : nx;
                size_t const actual_nu = out_degree(v, g) == 0 ? 0 : nu;
                size_t const actual_nc = out_degree(v, g) == 0 ? nct : nc;
                
                return OcpVertexSize {actual_nx, actual_nu, actual_nc, ns};
            }
        }
    {
    }


    DynamicOcpSize::DynamicOcpSize(size_t num_stages, size_t nx, size_t nu,
        size_t nc, size_t ns, size_t nct, bool first_state_empty)
    :	DynamicOcpSize {OcpTree(num_stages), nx, nu, nc, ns, nct, first_state_empty}
    {
    }
}