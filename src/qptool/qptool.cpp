#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include <nlohmann/json.hpp>

#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/core/Range.hpp>
#include <tmpc/core/PropertyMap.hpp>
//#include <tmpc/qp/HpmpcWorkspace.hpp>
//#include <tmpc/qp/TreeQpWorkspaceAdaptor.hpp>
#include <tmpc/qp/DualNewtonTreeWorkspace.hpp>
#include <tmpc/BlazeKernel.hpp>
#include <tmpc/json/JsonBlaze.hpp>
#include <tmpc/json/JsonQp.hpp>

#include <boost/graph/graphviz.hpp>


using namespace tmpc;


using Kernel = BlazeKernel<double>;
//using HpmpcSolver = TreeQpWorkspaceAdaptor<HpmpcWorkspace<Kernel>>;
using DualNewtonTreeSolver = DualNewtonTreeWorkspace<Kernel>;


template <typename QpSrc, typename QpDst>
inline void copyQpProperties(QpSrc const& src, QpDst& dst)
{
    auto const vert = vertices(src.graph());
    copyProperty(src.Q(), dst.Q(), vert);
    copyProperty(src.R(), dst.R(), vert);
    copyProperty(src.S(), dst.S(), vert);
    copyProperty(src.q(), dst.q(), vert);
    copyProperty(src.r(), dst.r(), vert);
    copyProperty(src.lx(), dst.lx(), vert);
    copyProperty(src.ux(), dst.ux(), vert);
    copyProperty(src.lu(), dst.lu(), vert);
    copyProperty(src.uu(), dst.uu(), vert);
    copyProperty(src.C(), dst.C(), vert);
    copyProperty(src.D(), dst.D(), vert);
    copyProperty(src.ld(), dst.ld(), vert);
    copyProperty(src.ud(), dst.ud(), vert);
    
    auto const edg = edges(src.graph());
    copyProperty(src.A(), dst.A(), edg);
    copyProperty(src.B(), dst.B(), edg);
    copyProperty(src.b(), dst.b(), edg);
}


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

    // Copy QP to the solver workspace.
    copyQpProperties(json_qp, solver);

    solver.print();

    // Solve the QP
    solver.solve();

    // Output the solution to json
    json j_sol;

    for (auto v : make_iterator_range(vertices(solver.graph())))
    {
        auto& node = j_sol["nodes"][v];
        node["x"] = get(solver.x(), v);
        node["u"] = get(solver.u(), v);

        json const input_node = j["nodes"][v];

        if (input_node.count("xopt"))
            node["x_opt"] = input_node["xopt"];

        if (input_node.count("uopt"))
            node["u_opt"] = input_node["uopt"];
    }

    // Print json to stdout
    std::cout << std::setw(4) << j_sol << std::endl;

    return EXIT_SUCCESS;
}