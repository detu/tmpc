#include <tmpc/integrator/ExplicitRungeKutta4.hpp>
#include <tmpc/integrator/MultiStepIntegrator.hpp>
#include <tmpc/test_tools.hpp>


namespace tmpc :: testing
{
	namespace 
	{
		blaze::StaticVector<double, 1> ode(double t, blaze::StaticVector<double, 1> const& x, blaze::StaticVector<double, 0> const& u)
		{
			return x;
		}
	}


	TEST(MultiStepIntegratorTest, testAccuracy)
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
}
