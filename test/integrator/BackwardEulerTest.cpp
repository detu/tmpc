#include <tmpc/integrator/ImplicitRungeKutta.hpp>
#include <tmpc/integrator/BackwardEuler.hpp>
#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	/// @brief Test integration of the problem dx/dt=x*u, x(0)=x_0 using backward Euler method
	///
	TEST(BackwardEulerTest, testIntegrateSimpleLinear)
	{
		using Real = double;
		size_t const NX = 1, NU = 1;
		using VecX = blaze::StaticVector<Real, NX, blaze::columnVector>;
		using VecU = blaze::StaticVector<Real, NU, blaze::columnVector>;
		
		ImplicitRungeKutta<Real> irk(NX, NU, backwardEuler<Real>());

		auto ode = [] (Real t, auto const& x, auto const& u, auto& f, auto& df_dx)
		{
			f = x * u;
			df_dx = {{u[0]}};
		};

		Real const t0 = 0.11;
		VecX const x0 {3.};
		VecU const u {2.};
		Real const h = 0.025;
		VecX const x1 = irk(ode, t0, x0, u, h);

		TMPC_EXPECT_APPROX_EQ(x1, x0 * exp(u * h), 0., 0.002);
	}


	/// @brief Test integration of a time-dependent ODE
	///
	TEST(BackwardEulerTest, testIntegrateTimeDependent)
	{
		using Real = double;
		size_t const NX = 1, NU = 1;
		using VecX = blaze::StaticVector<Real, NX, blaze::columnVector>;
		using VecU = blaze::StaticVector<Real, NU, blaze::columnVector>;
		
		ImplicitRungeKutta<Real> irk(NX, NU, backwardEuler<Real>());

		auto ode = [] (Real t, auto const& x, auto const& u, auto& f, auto& df_dx)
		{
			f = u * t;
			df_dx = {{0.}};
		};

		Real const t0 = 0.11;
		VecX const x0 {3.};
		VecU const u {2.};
		Real const h = 0.05;
		VecX const x1 = irk(ode, t0, x0, u, h);

		TMPC_EXPECT_APPROX_EQ(x1, x0 + u * (pow(t0 + h, 2) - pow(t0, 2)) / 2., 0., 0.001);
	}
}
