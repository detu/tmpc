#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include <nlohmann/json.hpp>

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/Range.hpp>
//#include <tmpc/qp/HpmpcWorkspace.hpp>
//#include <tmpc/qp/TreeQpWorkspaceAdaptor.hpp>
#include <tmpc/qp/DualNewtonTreeWorkspace.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/json/JsonBlaze.hpp>
#include <tmpc/json/JsonQp.hpp>
#include <tmpc/core/GraphTools.hpp>

#include <boost/graph/graphviz.hpp>


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

    using K = BlazeKernel<double>;
    JsonQp<K> json_qp(j);

    write_graphviz(std::cout, json_qp.graph());
    
    // Create solver workspace.
    DualNewtonTreeSolver solver {json_qp.graph(), json_qp.size()};

    solver.print();
    solver.solve();


    return EXIT_SUCCESS;
}