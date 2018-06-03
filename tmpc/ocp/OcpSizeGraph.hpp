#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/SizeT.hpp>

#include <boost/graph/adjacency_list.hpp>


namespace tmpc
{
    struct OcpVertex
    {
        OcpVertex() = default;


        OcpVertex(OcpSize const& sz)
        :   size(sz)
        {            
        }


        OcpSize size;
    };


    using OcpSizeGraph = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS, 
        OcpVertex,
        boost::property<boost::edge_index_t, size_t>
    >;
}
