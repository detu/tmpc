#include <tmpc/integrator/ExplicitRungeKutta4.hpp>
#include <tmpc/integrator/MultiStepIntegrator.hpp>

#include "PendulumData.hpp"

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	namespace 
	{
		blaze::StaticVector<double, 1> ode(double t, blaze::StaticVector<double, 1> const& x, blaze::StaticVector<double, 0> const& u)
		{
			return x;
		}
	}


	TEST(MultiStepIntegratorSimpleTest, testAccuracy)
	{
		double const h = 1.;
		double const h_max = 0.5;

		ExplicitRungeKutta4<double> rk4(1, 0);
		blaze::StaticVector<double, 1> const x0 {1.};
		blaze::StaticVector<double, 0> const u0 {};
		
		blaze::StaticVector<double, 1> const x1_true {exp(1.)};
		blaze::StaticVector<double, 1> const x1_rk4 = rk4(&ode, 0., x0, u0, h);

		MultiStepIntegrator<double> multistep(1, 0);
		blaze::StaticVector<double, 1> const x1_multistep = multistep(rk4, &ode, 0., x0, u0, h, h_max);

		double const err_rk4 = abs(x1_rk4[0] - x1_true[0]);
		double const err_multistep = abs(x1_multistep[0] - x1_true[0]);

		EXPECT_GE(err_rk4 / err_multistep, 0.6 * pow(h / h_max, 4));
	}


	class MultiStepIntegratorTest 
	: 	public Test
	{
	protected:
		using Real = double;
		typedef PendulumOdeBase ODE;

		PendulumOde ode_;
		PendulumOdeR odeR_;
		ExplicitRungeKutta4<Real> rk4_ {ODE::NX, ODE::NU};
		MultiStepIntegrator<Real> integrator_ {ODE::NX, ODE::NU};

		std::vector<TestPoint> test_data_;

		void SetUp() override
		{
			test_data_ = loadPendulumData();
			ASSERT_EQ(test_data_.size(), 600);
		}
	};


	TEST_F(MultiStepIntegratorTest, testIntegrate)
	{
		for (auto p : test_data_)
		{
			ODE::StateVector const xplus = integrator_(rk4_, ode_, p.t, p.x0, p.u, p.timeStep, p.timeStep * 0.1);

			// Tolerance increased by the factor of 1e+2, because the test data are obtained using 1-step RK4.
			ApproxEqual const is_approx(1e-8);
			EXPECT_PRED2(is_approx, forcePrint(xplus), forcePrint(p.xplus));
		}
	}


	TEST_F(MultiStepIntegratorTest, testIntegrateWithSensitivities)
	{
		for (auto p : test_data_)
		{
			ODE::StateVector xplus;
			ODE::StateStateMatrix A;
			ODE::StateInputMatrix B;
			integrator_(rk4_, ode_, p.t, p.x0, p.u, p.timeStep, p.timeStep * 0.1, xplus, A, B);

			// Tolerance increased by the factor of 1e+2, because the test data are obtained using 1-step RK4.
			EXPECT_PRED2(ApproxEqual(1e-8), forcePrint(xplus), forcePrint(p.xplus));
			EXPECT_PRED2(ApproxEqual(1e-8), forcePrint(A), forcePrint(p.A));
			EXPECT_PRED2(ApproxEqual(1e-8), forcePrint(B), forcePrint(p.B));
		}
	}
}
