/// @brief Demonstrates the use of the Newton solver

#include <tmpc/numeric/NewtonSolver.hpp>

#include <cmath>
#include <iostream>


/// @brief Functor for the residual function
struct Fun
{
    template <typename VT1, typename VT2, typename MT, bool SO>
    void operator()(blaze::Vector<VT1, blaze::columnVector> const& x, 
        blaze::Vector<VT2, blaze::columnVector>& f, blaze::Matrix<MT, SO>& J) const
    {
        ~f = {
            std::pow((~x)[0], 2) + std::pow((~x)[1], 3) - 1.,
            2. * (~x)[0] + 3. * std::pow((~x)[1], 2) - 4.
        };

        ~J = {
            {2. * (~x)[0], 3. * std::pow((~x)[1], 2)},
            {2., 6. * (~x)[1]}
        };
    };
};


/// @brief Iteration monitor
struct Monitor
{
    template <typename VT1, typename VT2, typename MT, bool SO>
    void operator()(size_t iter, 
        blaze::Vector<VT1, blaze::columnVector> const& x,
        blaze::Vector<VT2, blaze::columnVector> const& r,
        blaze::Matrix<MT, SO> const& J) const
    {
        std::cout << "iteration " << iter << std::endl;
        std::cout << "x = " << trans(~x);
        std::cout << "r = " << trans(~r);
        std::cout << "J = " << std::endl << J;
        std::cout << "----------------------------" << std::endl;
    }
};


int main(int, char **)
{
    size_t constexpr NX = 2;
        
    using Real = double;
    using Vec = blaze::StaticVector<Real, NX, blaze::columnVector>;
    using Mat = blaze::StaticMatrix<Real, NX, NX>;

    // Newton solver
    tmpc::NewtonSolver<Real> solver(NX);
    solver.maxIterations(20);

    // The function to minimize
    Fun fun;

    // Initial guess
    Vec const x0 {-2., -1.};
    
    // Find the solution
    Vec const x_star = solver.solve(fun, x0, Monitor());

    std::cout << "Solution found at " << trans(x_star);

    return 0;
}