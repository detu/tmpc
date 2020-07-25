#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/Exception.hpp>


namespace tmpc
{
    namespace 
    {
        auto scenarioTreeBranching(size_t num_stages, size_t robust_horizon, size_t branching)
        {
            using v_size_t = OcpTree::vertices_size_type;
            using d_size_t = OcpTree::degree_size_type;

            if (!(robust_horizon < num_stages))
                TMPC_THROW_EXCEPTION(std::invalid_argument("robust horizon must be less than number of stages"));
            
            std::vector<v_size_t> n_children;
            v_size_t current_stage_vertex_count = 1;

            for (size_t k = 0; k < num_stages; ++k)
            {
                v_size_t nc = k + 1 < num_stages ? (k < robust_horizon ? branching : 1) : 0;

                for (v_size_t i = 0; i < current_stage_vertex_count; ++i)
                    n_children.push_back(nc);

                current_stage_vertex_count *= nc;
            }

            return n_children;
        }        
    }


    OcpTree::OcpTree(size_t n_nodes)
    :   OcpTree(
            std::views::iota(
                OcpVertex {0},
                OcpVertex {static_cast<OcpTree::vertices_size_type>(n_nodes)}
            ) |
            std::views::transform(
                [n_nodes] (OcpVertex v) { return v + 1 < n_nodes ? 1 : 0; }
            )
        )
    {
    }


    OcpTree::OcpTree(size_t num_stages, size_t robust_horizon, size_t branching)
    :   OcpTree(scenarioTreeBranching(num_stages, robust_horizon, branching))
    {
    }
}