#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include <nlohmann/json.hpp>

#include <tmpc/ocp/OcpSizeGraph.hpp>
#include <tmpc/core/IteratorRange.hpp>
//#include <tmpc/qp/HpmpcWorkspace.hpp>
//#include <tmpc/qp/TreeQpWorkspaceAdaptor.hpp>
#include <tmpc/qp/DualNewtonTreeWorkspace.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/json/JsonBlaze.hpp>
#include <tmpc/json/JsonQp.hpp>
#include <tmpc/core/GraphTools.hpp>

#include <boost/range/adaptor/indexed.hpp>


using namespace tmpc;


using Kernel = BlazeKernel<double>;
//using HpmpcSolver = TreeQpWorkspaceAdaptor<HpmpcWorkspace<Kernel>>;
using DualNewtonTreeSolver = DualNewtonTreeWorkspace<Kernel>;


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

    for (auto const& j_vertex : j["nodes"] | indexed(0))
    {
        auto const& j_v = j_vertex.value();
        g[j_vertex.index()].size =
            OcpSize {j_v["q"].size(), j_v["r"].size(), j_v["ld"].size(), j_v["zl"].size()};
    }

    for (auto j_edge : j["edges"] | indexed(0))
        add_edge(j_edge.value()["from"], j_edge.value()["to"], j_edge.index(), g);


    // Create solver workspace.
    DualNewtonTreeSolver solver {g};
    //HpmpcSolver solver {g};

    // Set problem properties from json.
    /*
    for (auto const& j_vertex : j["nodes"] | indexed(0))
    {
        auto& qp_vertex = get(solver.problemVertex(), j_vertex.index());
        from_json(j_vertex.value(), qp_vertex);

        std::cout << qp_vertex.Q() << std::endl;
    }

    for (auto e : edgesR(solver.graph()))
    {
        auto const edge_index = get(solver.edgeIndex(), e);
        auto& qp_edge = get(solver.problemEdge(), e);
        from_json(j["edges"][edge_index], qp_edge);

        std::cout << qp_edge.A() << std::endl;
    }
    */

    DynamicMatrix<Kernel> Q {
        {1., 2., 3.}, 
        {3., 4., 5.},
        {6., 7., 8.}
    };
    
    DynamicMatrix<Kernel> R {
        {1., 2.}, {3., 4.}
    };

    DynamicMatrix<Kernel> S {
        {1., 2.}, 
        {3., 4.},
        {5., 6.}
    };

    DynamicVector<Kernel> q {11., 12., 13.};
    DynamicVector<Kernel> r {11., 12.};

    solver.setNodeObjective(0, Q, R, S, q, r);
    solver.print();
    solver.solve();


    return EXIT_SUCCESS;
}