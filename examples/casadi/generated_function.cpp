/// @brief Demonstrates the use of CasADi generated C functions from tmpc.
///
/// The sample function named "f" is defined in sample_functions.py

#include <tmpc/casadi/GeneratedFunction.hpp>
#include <tmpc/Math.hpp>

// Include the generated C header
#include <sample_functions.h>

#include <tuple>
#include <iostream>


int main(int, char **)
{
	// Create the GeneratedFunction object and bind it to the generated function named "f"
	tmpc::casadi::GeneratedFunction fun {f_functions()};

	// Function input arguments
	blaze::StaticMatrix<double, 3, 2> const A {
		{1., 2.},
		{3., 4.},
		{5., 6.}
	};
	blaze::StaticMatrix<double, 2, 2> const B {
		{7., 8.},
		{9., 10.}
	};	
	double const x = 0.1;

	// Function output arguments
	blaze::StaticMatrix<double, 3, 2> X;
	blaze::StaticVector<double, 2, blaze::rowVector> Y;

	// Call the CasADi generated function.
	// The operator() takes two tuples: the first tuple consists of input arguments,
	// the second tuple consists of output arguments.
	fun(std::tie(A, B, x), std::tie(X, Y));

	// Print the output and check it against the same formulas computed directly.
	std::cout << "X = \n" << X;
	std::cout << "(A * x) * B = \n" << evaluate((A * x) * B);
	std::cout << "Y = \n" << Y;
	std::cout << "[1, 1, 1] * (A * B) = \n" 
		<< evaluate(blaze::StaticVector<double, 3, blaze::rowVector> {1., 1., 1.} * (A * B));

	return 0;
}