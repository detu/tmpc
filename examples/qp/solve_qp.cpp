/// @brief Demonstrates how to create and solve an OCP QP.
///

#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/qp/Randomize.hpp>
#include <tmpc/qp/KktValue.hpp>
#include <tmpc/ocp/DynamicOcpSolution.hpp>
#include <tmpc/ocp/DynamicOcpKktValue.hpp>
#include <tmpc/print/qp/OcpQp.hpp>
#include <tmpc/print/ocp/OcpSolution.hpp>
#include <tmpc/print/ocp/OcpKktValue.hpp>
#include <tmpc/hpipm/NominalSolver.hpp>

#include <boost/exception/diagnostic_information.hpp>

#include <iostream>


int main(int, char **)
{
	using namespace tmpc;
	using Real = double;

	try
	{
		// Define and OCP QP graph with 3 vertices:
		//
		// (0)------>(1)------>(2)
		// v0   e0   v1   e1   v2
		//

		// Define the problem size at each node.
		// DynamicOcpSize {nx, nu} sets the state size to nx
		// and control size to nu at the corresponding node.
		DynamicOcpSize const size {
			{3, 2, 2},
			{2, 1, 1},
			{2, 0, 0}
		};

		// Construct the QP.
		// DynamicOcpQp represents an OCP QP
		// with size at each node defined at run-time.
		DynamicOcpQp<Real> qp {size};

		// Make a random feasible QP.
		randomize(qp);

		// Print the QP.
		std::cout << "**** QP ****" << std::endl;
		std::cout << qp;

		// Construct the solution object.
		// DynamicOcpSolution represents an OCP solution
		// with size at each node defined at run-time.
		DynamicOcpSolution<Real> sol {size};

		// Create the solver.
		// NominalSolver is the interface 
		// to the HPIPM solver https://github.com/giaf/hpipm
		// You need to compile and install HPIPM first.
		// Then set the CMake variable TMPC_WITH_HPIPM=ON when configuring tmpc.
		hpipm::NominalSolver<Real> solver {size};
		
		// Solve the QP.
		solver(qp, sol);

		// Print the solution.
		std::cout << "**** solution ****" << std::endl;
		std::cout << sol;

		// Check KKT values.
		DynamicOcpKktValue<Real> kkt {size};
		kktValue(qp, sol, kkt);

		std::cout << "**** KKT values ****" << std::endl;
		std::cout << kkt;
	}
	catch (boost::exception const& e)
	{
		// Catch errors and print diagnostic information, just in case.
		std::cerr << diagnostic_information(e);
	}

	return 0;
}