/// @brief Demonstrates the parallel evaluation of CasADi generated C functions from tmpc.
///
/// The sample function named "f" is defined in sample_functions.py

#include <tmpc/casadi/GeneratedFunction.hpp>
#include <tmpc/Math.hpp>

// Include the generated C header
#include <sample_functions.h>

#include <tuple>
#include <random>
#include <iostream>


int main(int, char **)
{
	// Create the GeneratedFunction object and bind it to the generated function named "f"
	tmpc::casadi::GeneratedFunction fun {f_functions()};
	
	// Number of function evaluations
	size_t const N = 5000;

	// Input arguments
	std::vector<blaze::StaticMatrix<double, 3, 2>> A(N);
	std::vector<blaze::StaticMatrix<double, 2, 2>> B(N);
	std::vector<double> x(N);

	// Randomize input arguments
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-1.0, 1.0);	

	for (size_t i = 0; i < N; ++i)
	{
		randomize(A[i]);
		randomize(B[i]);
		x[i] = dis(gen);
	}

	// Output arguments
	std::vector<blaze::StaticMatrix<double, 3, 2>> X(N), X_ref(N);
	std::vector<blaze::StaticVector<double, 2, blaze::rowVector>> Y(N), Y_ref(N);

	// Evaluate the function sequentially.
	for (size_t i = 0; i < N; ++i)
		fun(std::tie(A[i], B[i], x[i]), std::tie(X_ref[i], Y_ref[i]));

	// Evaluate the function in parallel.
	// Each thread gets its own copy of fun,
	// so the calculations in different threads do not interfere with each other!
	#pragma omp parallel num_threads(100) firstprivate(fun)
	{
		#pragma omp for
		for (size_t i = 0; i < N; ++i)
			fun(std::tie(A[i], B[i], x[i]), std::tie(X[i], Y[i]));
	}

	// Check if the results are the same.
	for (size_t i = 0; i < N; ++i)
		if (X[i] != X_ref[i] || Y[i] != Y_ref[i])
		{
			std::cerr << "mismatch at i=" << i << std::endl;
			return -1;
		}

	std::cout << "All OK!" << std::endl;

	return 0;
}