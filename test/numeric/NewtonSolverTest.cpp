#include <tmpc/numeric/NewtonSolver.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
    /// @brief Check that the Newton method finds the correct solution of a system of 2 equations.
    TEST(NewtonSolverTest, testSolve)
    {
        size_t constexpr NX = 2;
        
        using Real = double;
        using Vec = blaze::StaticVector<Real, NX, blaze::columnVector>;
        using Mat = blaze::StaticMatrix<Real, NX, NX>;

        // Newton solver
        NewtonSolver<Real> solver(NX);

        // Define the equation and its Jacobian
        auto fun = [] (auto const& x, auto& f, auto& J)
        {
            ~f = {
                pow(x[0], 2) + pow(x[1], 3) - 1.,
                2. * x[0] + 3. * pow(x[1], 2) - 4.
            };

            ~J = {
                {2. * x[0], 3. * pow(x[1], 2)},
                {2., 6. * x[1]}
            };
        };

        // Initial guess
        Vec const x0 {-2., -1.};

        // Find the solution
        Vec const x_star = solver.solve(fun, x0);
        
        // Check the solution
        TMPC_EXPECT_APPROX_EQ(x_star, (Vec {-2.48345, -1.72886}), 1.e-5, 0.);

        // Check residual at the solution
        Vec r;
        Mat J;
        fun(x_star, r, J);
        EXPECT_LT(maxNorm(r), solver.residualTolerance());
        EXPECT_EQ(solver.residualMaxNorm(), maxNorm(r));

        // Check number of iterations
        EXPECT_LE(solver.iterations(), solver.maxIterations());
    }


    /// @brief Solve the Rosenbrock problem https://en.wikipedia.org/wiki/Rosenbrock_function
    ///
    /// The corresponding system of equalities is
    /// 0 = -a + x + 2 b x^3 - 2 b x y
    /// 0 = b (-x^2 + y)
    ///
    /// The Jacobian of the system is
    /// 1 + 6 b x^2 - 2 b y, -2 b x
    /// -2 b x             , b
    ///
    /// Here we check how the line search works: if we did full Newton steps,
    /// then the redisuals would increase on iteration 2.
    /// So we check that the absolute value of the residual is strictly decreasing.
    TEST(NewtonSolverTest, testRosenbrock)
    {
        size_t constexpr NX = 2;
        
        using Real = double;
        using Vec = blaze::StaticVector<Real, NX, blaze::columnVector>;
        using Mat = blaze::StaticMatrix<Real, NX, NX>;

        // Newton solver
        NewtonSolver<Real> solver(NX);

        // Define the equation and its Jacobian
        auto fun = [] (auto const& x, auto& f, auto& J)
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

        // Initial guess
        Vec const x0 {-2., -1.};

        // Previous residual
        Vec r_prev(inf<Real>());

        // Find the solution
        Vec const x_star = solver.solve(fun, x0,
            [&r_prev] (size_t iter, auto const& x, auto const& r, auto const& J)
            {
                for (size_t i = 0; i < NX; ++i)
                    EXPECT_LT(abs(r[i]), abs(r_prev[i])) 
                        << " residual element " << i << " non-decreasing on iteration " << iter;
                r_prev = r;
            }
        );
        
        // Check the solution
        TMPC_EXPECT_APPROX_EQ(x_star, (Vec {1., 1.}), 1.e-10, 0.);
    }
}