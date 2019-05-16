#include <tmpc/integrator/ImplicitRungeKutta.hpp>
#include <tmpc/integrator/BackwardEuler.hpp>
#include <tmpc/integrator/GaussLegendre.hpp>
#include <tmpc/Testing.hpp>

#include <tuple>


namespace tmpc :: testing
{
	using Real = double;


	struct TestParam
	{
		ButcherTableau<Real> tableau;
		Real relTol;
		std::string methodName;
	};


	std::ostream& operator<<(std::ostream& os, TestParam const& p)
	{
		return os << p.methodName;
	}


	/// @brief Test integration of a simple linear ODE.
	///
	/// The first test param is the Butcher tableau, 
	/// the second is the relative tolerance for checking the integrator output.
	class IRK_SimpleLinearOdeTest
	:	public ::testing::TestWithParam<TestParam>
	{
	};


	/// @brief Test integration of a simple time-dependent ODE
	/// The first test param is the Butcher tableau, 
	/// the second is the relative tolerance for checking the integrator output.
	class IRK_SimpleTimeDependentOdeTest
	:	public ::testing::TestWithParam<TestParam>
	{
	};


	TEST_P(IRK_SimpleLinearOdeTest, testIntegrate)
	{
		using Real = double;
		size_t const NX = 1, NU = 1;
		using VecX = blaze::StaticVector<Real, NX, blaze::columnVector>;
		using VecU = blaze::StaticVector<Real, NU, blaze::columnVector>;
		
		ImplicitRungeKutta<Real> irk(NX, NU, GetParam().tableau);

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

		TMPC_EXPECT_APPROX_EQ(x1, x0 * exp(u * h), 0., GetParam().relTol);
	}

	
	TEST_P(IRK_SimpleTimeDependentOdeTest, testIntegrate)
	{
		using Real = double;
		size_t const NX = 1, NU = 1;
		using VecX = blaze::StaticVector<Real, NX, blaze::columnVector>;
		using VecU = blaze::StaticVector<Real, NU, blaze::columnVector>;
		
		ImplicitRungeKutta<Real> irk(NX, NU, GetParam().tableau);

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

		TMPC_EXPECT_APPROX_EQ(x1, x0 + u * (pow(t0 + h, 2) - pow(t0, 2)) / 2., 0., GetParam().relTol);
	}


	INSTANTIATE_TEST_SUITE_P(ImplicitRungeKuttaTest, IRK_SimpleLinearOdeTest,
		Values(
			TestParam {backwardEuler<Real>(), 0.002, "Backward Euler"},
			TestParam {gaussLegendre<Real>(2), 1e-7, "Gauss-Legendre"}
		)
	);


	INSTANTIATE_TEST_SUITE_P(ImplicitRungeKuttaTest, IRK_SimpleTimeDependentOdeTest,
		Values(
			TestParam {backwardEuler<Real>(), 0.001, "Backward Euler"},
			TestParam {gaussLegendre<Real>(2), 1e-7, "Gauss-Legendre"}
		)
	);
}
