#include <tmpc/ocp/OcpGraph.hpp>

#include <boost/iterator/iterator_facade.hpp>


namespace tmpc
{
    namespace
    {
        class LinearGraphEdgesIterator
        :   public boost::iterator_facade<
            LinearGraphEdgesIterator,
            std::pair<size_t, size_t> const,  // Value
            std::forward_iterator_tag  // CategoryOrTraversal
            //std::pair<size_t, size_t>   // Reference
        >
        {
        public:
            explicit LinearGraphEdgesIterator(size_t n)
            :   uv_(n, n + 1)
            {
            }


        private:
            friend class boost::iterator_core_access;


            std::pair<size_t, size_t> const& dereference() const
            {
                return uv_;
            }


            void increment() 
            { 
                ++uv_.first;
                ++uv_.second;
            }


            bool equal(LinearGraphEdgesIterator const& other) const
            {
                return uv_.first == other.uv_.first;
            }


            std::pair<size_t, size_t> uv_;
        };
    }


    OcpGraph ocpGraphLinear(size_t n_nodes)
    {
        // Create the graph.
        return n_nodes > 0 ? 
            OcpGraph(LinearGraphEdgesIterator(0), LinearGraphEdgesIterator(n_nodes - 1), n_nodes) : OcpGraph();
    }


    OcpGraph ocpGraphRobustMpc(size_t depth, size_t branching, size_t robust_horizon)
    {
        using v_size_t = OcpGraph::vertices_size_type;
        using d_size_t = OcpGraph::degree_size_type;
        
        std::vector<std::pair<v_size_t, v_size_t>> edge_list;
        std::vector<v_size_t> cur_layer, next_layer;

        v_size_t id_u = 0, id_v = 1;
        v_size_t target_vertex = 0;

        for (size_t k = 0; k < depth; ++k)
        {
            if (k == 0)
            {
                // Root node
                next_layer.push_back(target_vertex);
                ++target_vertex;
            }
            else
            {
                for (auto source_vertex : cur_layer)
                {
                    for (d_size_t j = 0; j < (k <= robust_horizon ? branching : 1); ++j)
                    {
                        next_layer.push_back(target_vertex);
                        edge_list.emplace_back(source_vertex, target_vertex);
                        ++target_vertex;
                    }
                }
            }

            swap(cur_layer, next_layer);
            next_layer.clear();
        }

        return OcpGraph(edge_list.begin(), edge_list.end(), target_vertex);
    }
}