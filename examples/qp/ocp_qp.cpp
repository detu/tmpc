/// @brief Demonstrates working with OCP QPs:
/// creating, accessing edge and vertex attributes, ...

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/qp/StaticOcpQp.hpp>
#include <tmpc/print/qp/OcpQp.hpp>

#include <iostream>
#include <array>


int main(int, char **)
{
    using namespace tmpc;
    using Real = double;

    // Number of states
    size_t constexpr NX = 3;

    // Number of controls
    size_t constexpr NU = 2;

    // Number of constraints
    size_t constexpr NC = 1;

    // Create an OCP tree from branching factors
    //
    //    v2---v4
    //   /
    // v0
    //   \
    //    v1---v3
    //
    OcpTree g {2,1,1,0,0};

    // Create a dynamically-sized OCP QP with specified
    // number of states, controls, and constraints for each vertex.
    DynamicOcpQp<Real> d_qp {
        DynamicOcpSize {g, 
            {
                {NX, NU, NC},
                {NX, NU, NC}, {NX, NU, NC},
                {NX, 0, NC}, {NX, 0, NC}
            }
        }
    };

    // Populate OCP QP properties
    //
    // The following line produces a compilation error, see
    // https://bitbucket.org/blaze-lib/blaze/issues/366
    // makePositiveDefinite(d_qp.Q(0));
    //
    d_qp.Q(0) = {
        {1., 2., 3.},
        {2., 4., 5.},
        {3., 5., 6.}};

    // Create a statically-sized OCP QP
    StaticOcpQp<Real, NX, NU, NC> s_qp {g};
    s_qp = d_qp;

    std::cout << s_qp << std::endl;

    return EXIT_SUCCESS;
}