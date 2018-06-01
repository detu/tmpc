#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/SizeT.hpp>

#include <boost/graph/adjacency_list.hpp>


namespace tmpc
{
    struct OcpSizeProperty_t
    {
        using kind = boost::vertex_property_tag;
    };


    using OcpSizeGraph = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::directedS, 
        boost::property<OcpSizeProperty_t, OcpSize>, 
        boost::property<boost::edge_index_t, size_t>
    >;
}
