#include <tmpc/integrator/ImplicitRungeKutta.hpp>
#include <tmpc/integrator/BackwardEulerMethod.hpp>
#include <tmpc/integrator/GaussLegendreMethod.hpp>

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	/// @brief Test integration of a simple time-dependent ODE.
	class IrkSimpleTimeDependentOdeTest
	:	public Test
	{
	protected:
		using Real = double;


		template <typename Method>
		void testIntegrate(Method const& method, double abs_tol, double rel_tol)
		{
			size_t const NX = 1, NZ = 0, NU = 1;
			using VecX = blaze::StaticVector<Real, NX, blaze::columnVector>;
			using VecU = blaze::StaticVector<Real, NU, blaze::columnVector>;
			
			ImplicitRungeKutta<Real> irk(method, NX, NZ, NU);

			auto ode = [] (Real t, auto const& xdot, auto const& x, auto const& z,
				auto const& u, auto& f, auto& Jxdot, auto& Jx, auto& Jz)
			{
				f = u * t - xdot;
				Jxdot = {{-1.}};
				Jx = {{0.}};
			};

			Real const t0 = 0.11;
			VecX const x0 {3.};
			VecU const u {2.};
			Real const h = 0.05;
			VecX x1;
			irk(ode, t0, h, x0, u, x1);

			TMPC_EXPECT_APPROX_EQ(x1, x0 + u * (pow(t0 + h, 2) - pow(t0, 2)) / 2., abs_tol, rel_tol);
		}
	};
	

	TEST_F(IrkSimpleTimeDependentOdeTest, testBackwardEuler)
	{
		testIntegrate(BackwardEulerMethod(), 0., 0.001);
	}


	TEST_F(IrkSimpleTimeDependentOdeTest, testGaussLegendre2)
	{
		testIntegrate(GaussLegendreMethod(2), 0., 1e-14);
	}


	TEST_F(IrkSimpleTimeDependentOdeTest, testGaussLegendre3)
	{
		testIntegrate(GaussLegendreMethod(3), 0., 1e-17);
	}
}