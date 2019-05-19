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
    solver.maxIterations(1000);
    solver.backTrackingAlpha(0.7);

    // Define the equation and its Jacobian.
    // This is the Rosenbrock problem.
    auto rosenbrock = [] (auto const& x, auto& f, auto& J)
    {
        Real a = 1.;
        Real b = 10.;

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
    Vec const x_star = solver.solve(rosenbrock, x0, monitor);
    std::cout << "Rosenbrock problem solution: " << trans(x_star);

    return 0;
}