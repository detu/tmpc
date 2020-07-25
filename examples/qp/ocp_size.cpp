/// @brief Demonstrates working with OCP size objects.

#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/ocp/StaticOcpSize.hpp>
#include <tmpc/print/ocp/OcpSize.hpp>

#include <iostream>


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
    OcpTree const g {2,1,1,0,0};

    // Create a dynamically-sized OCP size with specified
    // number of states, controls, and constraints for each vertex.
    DynamicOcpSize size1 {g, 
        {
            {NX, NU, NC},
            {NX, NU, NC}, {NX, NU, NC},
            {NX, 0, NC}, {NX, 0, NC}
        }
    };

    // Create an equivalent statically-sized OCP size.
    StaticOcpSize<NX, NU, NC> size2 {g};

    // Print sizes.
    std::cout << "size1:" << std::endl << size1;
    std::cout << "size2:" << std::endl << size2;

    // Compare sizes.
    std::cout << "Sizes are equal: " << (size1 == size2 ? "yes" : "no") << std::endl;

    return EXIT_SUCCESS;
}