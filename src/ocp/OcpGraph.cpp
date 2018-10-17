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
}