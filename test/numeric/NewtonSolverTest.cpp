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
        Vec x_star;
        solver(fun, x0, x_star);
        
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


    /// @brief Test solution sensitivities correctness for scalar equation 0=x^2-p
    TEST(NewtonSolverTest, testSolutionSensitivityScalarQuadratic)
    {
        size_t constexpr NX = 1;
        size_t constexpr NP = 1;

        using Real = double;
        using VecX = blaze::StaticVector<Real, NX, blaze::columnVector>;
        using VecP = blaze::StaticVector<Real, NP, blaze::columnVector>;
        using MatXX = blaze::StaticMatrix<Real, NX, NX, blaze::columnMajor>;
        using MatXP = blaze::StaticMatrix<Real, NX, NP, blaze::columnMajor>;
        Real const p = 2.;
        
        // Define the equation and its Jacobian
        auto fun = [p] (auto const& x, auto& f, auto& J)
        {
            ~f = {pow(x[0], 2) - p};
            ~J = {{2. * x[0]}};
        };

        // The derivative of the equation w.r.t. the parameter
        auto dfdp = [p] (auto const& x, auto& df_dp)
        {
            ~df_dp = {{-1.}};
        };

        // Newton solver
        NewtonSolver<Real> solver(NX);

        // Initial guess
        VecX const x0 {1.};

        // Find the solution and its sensitivities
        VecX x_star;
        MatXX x_sens;
        solver(fun, dfdp, x0, x_star, x_sens);
        
        // Check the solution
        TMPC_EXPECT_APPROX_EQ(x_star, (VecX {sqrt(p)}), 1.e-10, 0.);

        // Check residual at the solution
        VecX r;
        MatXX J;
        fun(x_star, r, J);
        EXPECT_LT(maxNorm(r), solver.residualTolerance());
        EXPECT_EQ(solver.residualMaxNorm(), maxNorm(r));

        // Check solution sensitivity
        TMPC_EXPECT_APPROX_EQ(x_sens, (MatXP {{1. / (2. * sqrt(p))}}), 1.e-10, 0.);

        // Check number of iterations
        EXPECT_LE(solver.iterations(), solver.maxIterations());
    }


    /// @brief Test solution sensitivities correctness for vector equation 0=A*x-p
    TEST(NewtonSolverTest, testSolutionSensitivityVectorLinear)
    {
        size_t constexpr NX = 2;
        size_t constexpr NP = 2;

        using Real = double;
        using VecX = blaze::StaticVector<Real, NX, blaze::columnVector>;
        using VecP = blaze::StaticVector<Real, NP, blaze::columnVector>;
        using MatXX = blaze::StaticMatrix<Real, NX, NX, blaze::columnMajor>;
        using MatXP = blaze::StaticMatrix<Real, NX, NP, blaze::columnMajor>;
        
        MatXX const A {
            {1., 2.},
            {-3., 4.}
        };

        VecP const p {1., 0.5};
        
        // Define the equation and its Jacobian
        auto fun = [&A, &p] (auto const& x, auto& f, auto& J)
        {
            ~f = A * ~x - p;
            ~J = A;
        };

        // The derivative of the equation w.r.t. the parameter
        auto dfdp = [p] (auto const& x, auto& df_dp)
        {
            ~df_dp = -blaze::IdentityMatrix<Real>(NX);
        };

        // Newton solver
        NewtonSolver<Real> solver(NX);

        // Initial guess
        VecX const x0 {0., 0.};

        // Find the solution and its sensitivities
        VecX x_star;
        MatXX x_sens;
        solver(fun, dfdp, x0, x_star, x_sens);
        
        // Check the solution
        TMPC_EXPECT_APPROX_EQ(x_star, evaluate(inv(A) * p), 1.e-10, 0.);

        // Check residual at the solution
        VecX r;
        MatXX J;
        fun(x_star, r, J);
        EXPECT_LT(maxNorm(r), solver.residualTolerance());
        EXPECT_EQ(solver.residualMaxNorm(), maxNorm(r));

        // Check solution sensitivity
        TMPC_EXPECT_APPROX_EQ(x_sens, evaluate(inv(A)), 1.e-10, 0.);

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
    ///
    /// NOTE:
    /// The test reproduces the issue https://gitlab.syscop.de/mikhail.katliar/tmpc/issues/52
    /// and is disabled until the issue is fixed.
    /// 
    TEST(NewtonSolverTest, DISABLED_testRosenbrockLineSearch)
    {
        size_t constexpr NX = 2;
        
        using Real = double;
        using Vec = blaze::StaticVector<Real, NX, blaze::columnVector>;
        using Mat = blaze::StaticMatrix<Real, NX, NX>;

        // Newton solver
        NewtonSolver<Real> solver(NX);
        solver.maxIterations(20000);
        solver.backtrackingAlpha(0.5);

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
        Vec x_star;
        solver(fun, x0, x_star,
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
	
	
    /// @brief Solve the Rosenbrock problem https://en.wikipedia.org/wiki/Rosenbrock_function
    /// for multiple random initial points.
    ///
    /// The corresponding system of equalities is
    /// 0 = -a + x + 2 b x^3 - 2 b x y
    /// 0 = b (-x^2 + y)
    ///
    /// The Jacobian of the system is
    /// 1 + 6 b x^2 - 2 b y, -2 b x
    /// -2 b x             , b
    ///
    /// Here we test that the algorithm can find the solution starting from different initial points.
    TEST(NewtonSolverTest, testRosenbrockMultipleInitialPoints)
    {
        size_t constexpr NX = 2;
        
        using Real = double;
        using Vec = blaze::StaticVector<Real, NX, blaze::columnVector>;
        using Mat = blaze::StaticMatrix<Real, NX, NX>;

        // Newton solver
        NewtonSolver<Real> solver(NX);
        solver.residualTolerance(1e-10);
        solver.maxIterations(20);

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

        // True solution
        Vec const x_true {1., 1.};

        // Try different initial points
        size_t const num_points = 10000;
        for (size_t count = 0; count < num_points; ++count)
        {
            // Initial guess
            Vec x0;
            randomize(x0, -20., 20.);

            // Find the solution
            Vec x_star;
            ASSERT_NO_THROW(solver(fun, x0, x_star))
                << " at starting point " << trans(x0);
            
            // Check the solution
            TMPC_EXPECT_APPROX_EQ(x_star, x_true, 1e-9, 0.)
                << " the difference is " << trans(x_star - x_true);
        }
    }
}