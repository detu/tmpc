/// @brief Demonstrates how to convert an OCP QP to JSON.
///

#include <tmpc/json/JsonQp.hpp>
#include <tmpc/json/JsonBlaze.hpp>

#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <iostream>


int main(int, char **)
{
    using namespace tmpc;

    using Kernel = BlazeKernel<double>;

    JsonQp<Kernel> qp;
    auto const v0 = add_vertex(OcpSize {3, 2, 0}, qp.graph());
    auto const v1 = add_vertex(OcpSize {2, 1, 0}, qp.graph());
    auto const e0 = add_edge(v0, v1, 0, qp.graph()).first;

    put(qp.Q(), v0, DynamicMatrix<Kernel> {3, 3});
    put(qp.q(), v0, DynamicVector<Kernel>(3u, 0.));
    put(qp.A(), e0, DynamicMatrix<Kernel> {2, 3});

    std::cout << std::setw(4) << qp.json() << std::endl;

    std::cout << get(qp.Q(), v0) << std::endl;
    std::cout << get(qp.q(), v0) << std::endl;

    /*
    using Stage = OcpQp<Kernel>;

    Stage stage0 {OcpSize {3, 2, 0}, 2};
    Stage stage1 {OcpSize {2, 1, 0}, 0};
    
    randomize(stage0);
    randomize(stage1);

    std::vector<Stage> qp;
    qp.push_back(stage0);
    qp.push_back(stage1);

    json j = ocpQpToJson(qp);

    std::cout << std::setw(4) << j << std::endl;
    */
    
    return 0;
}
