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
    solver.maxIterations(10);

    // Define the equation and its jacobian
    auto fun = [] (auto const& x, auto& f, auto& J)        
    {
        ~f = {
            std::pow(x[0], 2) + std::pow(x[1], 3) - 1.,
            2. * x[0] + 3. * std::pow(x[1], 2) - 4.
        };

        ~J = {
            {2. * x[0], 3. * std::pow(x[1], 2)},
            {2., 6. * x[1]}
        };
    };

    // Initial point
    Vec const x0 {-2., -1.};

     // Monitor function
    auto monitor = [] (size_t iter, auto const& x, auto const& r, auto const& J)
    {
        std::cout << "iteration " << iter << std::endl;
        std::cout << "x = " << trans(~x);
        std::cout << "r = " << trans(~r);
        std::cout << "J = " << std::endl << J;
        std::cout << "----------------------------" << std::endl;
    };

    // Find the solution
    Vec const x_star = solver.solve(fun, x0, monitor);
    std::cout << "Solution found at " << trans(x_star);

    return 0;
}