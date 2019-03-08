#include <tmpc/integrator/ExplicitRungeKutta4.hpp>

#include "PendulumOde.hpp"
#include "PendulumData.hpp"

#include <tmpc/test_tools.hpp>

#include <fstream>


namespace tmpc :: testing
{
	class ExplicitRungeKutta4Test 
	: 	public ::testing::Test
	{
	protected:
		using Real = double;
		typedef PendulumOdeBase ODE;
		typedef ExplicitRungeKutta4<Real> Integrator;

		PendulumOde ode_;
		PendulumOdeR odeR_;
		Integrator integrator_ {ODE::NX, ODE::NU};

		std::vector<TestPoint> test_data_;

		void SetUp() override
		{
			test_data_ = loadPendulumData();
			ASSERT_EQ(test_data_.size(), 600);
		}
	};


	TEST_F(ExplicitRungeKutta4Test, testIntegrate)
	{
		for (auto p : test_data_)
		{
			ODE::StateVector const xplus = integrator_(ode_, p.t, p.x0, p.u, p.timeStep);

			MatrixApproxEquality const is_approx(1e-10);
			EXPECT_PRED2(is_approx, xplus, p.xplus);
		}
	}


	TEST_F(ExplicitRungeKutta4Test, testIntegrateWithSensitivities)
	{
		for (auto p : test_data_)
		{
			ODE::StateVector xplus;
			ODE::StateStateMatrix A;
			ODE::StateInputMatrix B;
			integrator_(ode_, p.t, p.x0, p.u, p.timeStep, xplus, A, B);

			/*
			EXPECT_EQ(forcePrint(xplus), forcePrint(p.xplus));
			EXPECT_EQ(forcePrint(A), forcePrint(p.A));
			EXPECT_EQ(forcePrint(B), forcePrint(p.B));
			*/

			MatrixApproxEquality const is_approx(1e-10);
			EXPECT_PRED2(is_approx, xplus, p.xplus);
			EXPECT_PRED2(is_approx, A, p.A);
			EXPECT_PRED2(is_approx, B, p.B);
		}
	}
}


// TEST_F(ExplicitRungeKutta4Test, integrate_q_correct)
// {
// 	TestPoint p;

// 	unsigned count = 0;
// 	while (test_data_ >> p)
// 	{
// 		ODE::StateVector xplus;
// 		ODE::StateStateMatrix A;
// 		ODE::StateInputMatrix B;
// 		ODE::QuadVector qf;
// 		ODE::QuadStateMatrix qA;
// 		ODE::QuadInputMatrix qB;
// 		integrate(integrator_, ode_, p.t, p.x0, p.u, xplus, A, B, qf, qA, qB);

// 		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));
// 		EXPECT_THAT(as_container(A    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.A    )));
// 		EXPECT_THAT(as_container(B    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.B    )));
// 		EXPECT_THAT(as_container(qf   ), testing::Pointwise(FloatNearPointwise(1e-4), as_container(p.qf   )));
// 		EXPECT_THAT(as_container(qA   ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.qA   )));
// 		EXPECT_THAT(as_container(qB   ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.qB   )));

// 		++count;
// 	}

// 	EXPECT_EQ(count, 600);
// }

// TEST_F(ExplicitRungeKutta4Test, integrate_r_correct)
// {
// 	TestPoint p;

// 	unsigned count = 0;
// 	while (test_data_ >> p)
// 	{
// 		ODE::StateVector xplus;
// 		ODE::StateStateMatrix A;
// 		ODE::StateInputMatrix B;
// 		double cf;
// 		ODE::StateVector cA;
// 		ODE::InputVector cB;
// 		ODE::StateStateMatrix cQ;
// 		ODE::InputInputMatrix cR;
// 		ODE::StateInputMatrix cS;
// 		integrate(integrator_, odeR_, p.t, p.x0, p.u, xplus, A, B, cf, cA, cB, cQ, cR, cS);

// 		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));
// 		EXPECT_THAT(as_container(A    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.A    )));
// 		EXPECT_THAT(as_container(B    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.B    )));
// 		EXPECT_NEAR(cf, p.cf, 1e-10);
// 		EXPECT_THAT(as_container(cA   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cA   )));
// 		EXPECT_THAT(as_container(cB   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cB   )));
// 		EXPECT_THAT(as_container(cQ   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cQ   )));
// 		EXPECT_THAT(as_container(cR   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cR   )));
// 		EXPECT_THAT(as_container(cS   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cS   )));
// 		++count;
// 	}

// 	EXPECT_EQ(count, 600);
// }

// TEST_F(ExplicitRungeKutta4Test, integrate_qr_correct)
// {
// 	TestPoint p;

// 	unsigned count = 0;
// 	while (test_data_ >> p)
// 	{
// 		ODE::StateVector xplus;
// 		ODE::StateStateMatrix A;
// 		ODE::StateInputMatrix B;
// 		ODE::QuadVector qf;
// 		ODE::QuadStateMatrix qA;
// 		ODE::QuadInputMatrix qB;
// 		double cf;
// 		ODE::StateVector cA;
// 		ODE::InputVector cB;
// 		ODE::StateStateMatrix cQ;
// 		ODE::InputInputMatrix cR;
// 		ODE::StateInputMatrix cS;
// 		integrate(integrator_, ode_, p.t, p.x0, p.u, xplus, A, B, qf, qA, qB, cf, cA, cB, cQ, cR, cS);

// 		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));
// 		EXPECT_THAT(as_container(A    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.A    )));
// 		EXPECT_THAT(as_container(B    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.B    )));
// 		EXPECT_THAT(as_container(qf   ), testing::Pointwise(FloatNearPointwise(1e-4), as_container(p.qf   )));
// 		EXPECT_THAT(as_container(qA   ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.qA   )));
// 		EXPECT_THAT(as_container(qB   ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.qB   )));
// 		EXPECT_NEAR(cf, p.cf, 1e-10);
// 		EXPECT_THAT(as_container(cA   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cA   )));
// 		EXPECT_THAT(as_container(cB   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cB   )));
// 		EXPECT_THAT(as_container(cQ   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cQ   )));
// 		EXPECT_THAT(as_container(cR   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cR   )));
// 		EXPECT_THAT(as_container(cS   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cS   )));
// 		++count;
// 	}

// 	EXPECT_EQ(count, 600);
// }

// TEST_F(ExplicitRungeKutta4Test, integrate_no_sens_correct)
// {
// 	TestPoint p;

// 	unsigned count = 0;
// 	while (test_data_ >> p)
// 	{
// 		auto const xplus = integrate(integrator_, ode_, p.t, p.x0, p.u);
// 		EXPECT_PRED2(MatrixApproxEquality(1e-5), xplus, p.xplus);

// 		++count;
// 	}

// 	EXPECT_EQ(count, 600);
// }

// TEST_F(ExplicitRungeKutta4Test, integrate_fd_correct)
// {
// 	TestPoint p;

// 	unsigned count = 0;
// 	while (test_data_ >> p)
// 	{
// 		ODE::StateVector xplus;
// 		ODE::StateStateMatrix A;
// 		ODE::StateInputMatrix B;
// 		integrate(integrator_, ode_, p.t, p.x0, p.u, xplus, A, B);

// 		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));
// 		EXPECT_THAT(as_container(A), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.A)));
// 		EXPECT_THAT(as_container(B), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.B)));

// 		++count;
// 	}

// 	EXPECT_EQ(count, 600);
// }
