#include <tmpc/numeric/NewtonSolver.hpp>
#include <tmpc/Testing.hpp>

#include <cmath>


namespace tmpc :: testing
{
    class NewtonSolverTest
    :   public Test
    {
    protected:
        // Test functor
        struct Fun
        {
            template <typename VT1, typename VT2, typename MT, bool SO>
            void operator()(blaze::Vector<VT1, blaze::columnVector> const& x, 
                blaze::Vector<VT2, blaze::columnVector>& f, blaze::Matrix<MT, SO>& J) const
            {
                ~f = {
                    pow((~x)[0], 2) + pow((~x)[1], 3) - 1.,
                    2. * (~x)[0] + 3. * pow((~x)[1], 2) - 4.
                };

                ~J = {
                    {2. * (~x)[0], 3. * pow((~x)[1], 2)},
                    {2., 6. * (~x)[1]}
                };
            };
        };


        Fun fun_;
    };


    /// @brief Check that the Newton method finds the correct solution of a system of 2 equations.
    TEST_F(NewtonSolverTest, testSolve)
    {
        size_t constexpr NX = 2;
        
        using Real = double;
        using Vec = blaze::StaticVector<Real, NX, blaze::columnVector>;
        using Mat = blaze::StaticMatrix<Real, NX, NX>;

        // Newton solver
        NewtonSolver<Real> solver(NX);
        solver.maxIterations(20);

        // Initial guess
        Vec const x0 {-2., -1.};

        // Find the solution
        Vec const x_star = solver.solve(fun_, x0);
        
        // Check the solution
        TMPC_EXPECT_APPROX_EQ(x_star, (Vec {-2.48345, -1.72886}), 1.e-5, 0.);

        // Check residual at the solution
        Vec r;
        Mat J;
        fun_(x_star, r, J);
        EXPECT_LT(maxNorm(r), solver.residualTolerance());
        EXPECT_EQ(solver.residualMaxNorm(), maxNorm(r));

        // Check number of iterations
        EXPECT_LE(solver.iterations(), solver.maxIterations());
    }
}