/// @brief Demonstrates how to convert an OCP QP to JSON.
///

#include <tmpc/json/JsonQp.hpp>
#include <tmpc/json/JsonBlaze.hpp>

#include <tmpc/qp/OcpQpStage.hpp>
#include <tmpc/BlazeKernel.hpp>

#include <iostream>


int main(int, char **)
{
    using namespace tmpc;

    using Kernel = BlazeKernel<double>;

    std::vector<DynamicOcpSize> size =
    {
        DynamicOcpSize {3, 2, 0},
        DynamicOcpSize {2, 1, 0}
    };

    JsonQp<Kernel> qp;
    auto const v0 = add_vertex(qp.graph());
    auto const v1 = add_vertex(qp.graph());
    auto const e0 = add_edge(v0, v1, qp.graph()).first;

    put(qp.Q(), v0, DynamicMatrix<Kernel> {3, 3});
    put(qp.q(), v0, DynamicVector<Kernel>(3u, 0.));
    put(qp.A(), e0, DynamicMatrix<Kernel> {2, 3});

    std::cout << std::setw(4) << qp.json() << std::endl;

    std::cout << get(qp.Q(), v0) << std::endl;
    std::cout << get(qp.q(), v0) << std::endl;

    /*
    using Stage = OcpQpStage<Kernel>;

    Stage stage0 {DynamicOcpSize {3, 2, 0}, 2};
    Stage stage1 {DynamicOcpSize {2, 1, 0}, 0};
    
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
