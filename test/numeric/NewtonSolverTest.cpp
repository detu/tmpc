#include <tmpc/numeric/NewtonSolver.hpp>
#include <tmpc/Testing.hpp>

#include <cmath>


namespace tmpc :: testing
{
    using std::cos;
    using std::sin;
    using std::pow;
    using std::exp;


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
                    {2. * (~x)[0], 3. * (~x)[1]},
                    {2., 6. * (~x)[1]}
                };
            };
        };


        Fun fun_;
    };


    /// @brief Check that the Newton method finds the correct solution of a system of 2 equations.
    /// 
    /// The example is taken from here: https://www.mathworks.com/help/optim/ug/fsolve.html
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
        Vec const x0 {-2.48345, -1.72886}; //{-2., -1.};
        
        // You will be surprised, by they are really equal within machine precision.
        TMPC_EXPECT_APPROX_EQ(solver.solve(fun_, x0), (Vec {-2.48345, -1.72886}), 1.e-5, 0.);
    }
}