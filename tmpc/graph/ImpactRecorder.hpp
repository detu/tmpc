#pragma once

#include "GraphTraits.hpp"
#include "Visitors.hpp"
#include "DepthFirstSearch.hpp"


namespace tmpc :: graph
{
    template <typename ImpactMap>
    class ImpactRecorder
    :   public base_visitor<ImpactRecorder<ImpactMap>>
    {
    public:
        using event_filter = on_finish_vertex;


        ImpactRecorder(ImpactMap impact)
        :   impact_(impact)
        {
        }


        template <class Vertex, class Graph>
        void operator()(Vertex u, const Graph& g) 
        {
            auto const children = adjacent_vertices(u, g);

            if (children.empty())
                put(impact_, u, 1);
            else
            {
                typename graph_traits<Graph>::vertices_size_type n = 0;
                for (auto v : children)
                    n += get(impact_, v);

                put(impact_, u, n);
            }
        }
        

    private:
        ImpactMap impact_;
    };
}