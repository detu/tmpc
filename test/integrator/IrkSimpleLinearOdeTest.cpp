#include <tmpc/integrator/ImplicitRungeKutta.hpp>
#include <tmpc/integrator/BackwardEulerMethod.hpp>
#include <tmpc/integrator/GaussLegendreMethod.hpp>

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	/// @brief Test integration of a simple linear ODE.
	class IrkSimpleLinearOdeTest
	:	public Test
	{
	protected:
		template <typename Method>
		void testIntegrate_impl(Method const& method, double abs_tol, double rel_tol)
		{
			using Real = double;
			size_t const NX = 1, NZ = 0, NU = 1;
			using VecX = blaze::StaticVector<Real, NX, blaze::columnVector>;
			using VecU = blaze::StaticVector<Real, NU, blaze::columnVector>;
			
			ImplicitRungeKutta<Real> irk(method, NX, NZ, NU);

			auto ode = [] (Real t, auto const& xdot, auto const& x, auto const& z,
				auto const& u, auto& f, auto& Jxdot, auto& Jx, auto& Jz)
			{
				f = x * u - xdot;
				Jxdot = {{-1.}};
				Jx = {{u[0]}};
			};

			Real const t0 = 0.11;
			VecX const x0 {3.};
			VecU const u {2.};
			Real const h = 0.025;
			VecX x1;
			irk(ode, t0, h, x0, u, x1);

			TMPC_EXPECT_APPROX_EQ(x1, x0 * exp(u * h), abs_tol, rel_tol);
		}
	};


	TEST_F(IrkSimpleLinearOdeTest, testBackwardEuler)
	{
		testIntegrate_impl(BackwardEulerMethod(), 0., 0.002);
	}


	TEST_F(IrkSimpleLinearOdeTest, testGaussLegendre2)
	{
		testIntegrate_impl(GaussLegendreMethod(2), 0., 1e-7);
	}


	TEST_F(IrkSimpleLinearOdeTest, testGaussLegendre3)
	{
		testIntegrate_impl(GaussLegendreMethod(3), 0., 1e-14);
	}
}