#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include <nlohmann/json.hpp>

#include <tmpc/ocp/OcpSizeGraph.hpp>
#include <tmpc/core/IteratorRange.hpp>
#include <tmpc/qp/HpmpcWorkspace.hpp>
#include <tmpc/qp/TreeQpWorkspaceAdaptor.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/json/JsonBlaze.hpp>
#include <tmpc/json/JsonQp.hpp>

#include <boost/range/adaptor/indexed.hpp>


using namespace tmpc;


using Kernel = BlazeKernel<double>;
using HpmpcSolver = TreeQpWorkspaceAdaptor<HpmpcWorkspace<Kernel>>;


int main(int argc, char ** argv)
{
    using nlohmann::json;
    using boost::adaptors::indexed;

    json j;
    if (argc > 1)
        std::ifstream(argv[1]) >> j;
    else
        std::cin >> j;

    // Build size graph from json.
    size_t const N = j["nodes"].size();
    OcpSizeGraph g(N);

    auto ocp_size = get(OcpSizeProperty_t(), g);
    for (auto const& j_vertex : j["nodes"] | indexed(0))
    {
        auto const& j_v = j_vertex.value();
        put(ocp_size, j_vertex.index(), 
            OcpSize {j_v["q"].size(), j_v["r"].size(), j_v["ld"].size(), j_v["zl"].size()});
    }

    for (auto j_edge : j["edges"] | indexed(0))
        add_edge(j_edge.value()["from"], j_edge.value()["to"], j_edge.index(), g);


    // Create solver workspace.
    HpmpcSolver solver {g};

    // Set problem properties from json.
    for (auto const& j_vertex : j["nodes"] | indexed(0))
    {
        auto& qp_vertex = get(solver.problemVertex(), j_vertex.index());
        from_json(j_vertex.value(), qp_vertex);

        std::cout << qp_vertex.Q() << std::endl;
    }


    return EXIT_SUCCESS;
}