#pragma once

#include <boost/graph/visitors.hpp>


namespace tmpc :: graph
{
    using boost::on_finish_vertex;
    using boost::on_discover_vertex;
    using boost::on_tree_edge;

    using boost::base_visitor;
    using boost::distance_recorder;
    using boost::record_distances;
}