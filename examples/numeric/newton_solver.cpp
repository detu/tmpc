/// @brief Demonstrates the use of the Newton solver

#include <tmpc/numeric/NewtonSolver.hpp>

#include <iostream>


int main(int, char **)
{
    size_t constexpr NX = 2;
        
    using Real = double;
    using Vec = blaze::StaticVector<Real, NX, blaze::columnVector>;
    using Mat = blaze::StaticMatrix<Real, NX, NX>;

    // Newton solver
    tmpc::NewtonSolver<Real> solver(NX);
    solver.maxIterations(200000);
    solver.backtrackingAlpha(0.5);

    // Define the equation and its Jacobian.
    // This is the Rosenbrock problem.
    auto rosenbrock = [] (auto const& x, auto& f, auto& J)
    {
        Real a = 1.;
        Real b = 100.;

        f = {
            -a + x[0] + 2. * b * pow(x[0], 3) - 2. * b * x[0] * x[1],
            b * (-pow(x[0], 2) + x[1])
        };

        J = {
            {1. + 6. * b * pow(x[0], 2) - 2. * b * x[1], -2. * b * x[0]},
            {-2. * b * x[0], b}
        };
    };

    // Initial point
    Vec const x0 {-2., -1.};

     // Monitor function
    auto monitor = [] (size_t iter, auto const& x, auto const& r, auto const& J)
    {
        std::cout << "iteration " << iter << std::endl;
        std::cout << "x = " << trans(x);
        std::cout << "r = " << trans(r);
        std::cout << "J = " << std::endl << J;
        std::cout << "----------------------------" << std::endl;
    };

    // Find the solution
    Vec x_star;
    solver(rosenbrock, x0, x_star, monitor);
    std::cout << "Rosenbrock problem solution: " << trans(x_star);
    std::cout << "Total number of Newton iterations: " << solver.iterations() << std::endl;
    std::cout << "Total number of function evaluations: " << solver.functionEvaluations() << std::endl;
    std::cout << "Average number of function evaluations per Newton iteration: " 
        << static_cast<double>(solver.functionEvaluations()) / solver.iterations() << std::endl;

    return 0;
}