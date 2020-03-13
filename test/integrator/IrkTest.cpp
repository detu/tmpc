#include <tmpc/integrator/ImplicitRungeKutta.hpp>
#include <tmpc/integrator/BackwardEulerMethod.hpp>

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	TEST(IrkTest, testInvalidDimensionX)
	{
		using Real = double;

		size_t const NX = 1, NZ = 0, NU = 1;
		ImplicitRungeKutta<Real> irk(BackwardEulerMethod(), NX, NZ, NU);

		auto foo_ode = [] (Real t, auto const& xdot, auto const& x, auto const& z, 
			auto const& u, auto& f, auto& Jxdot, auto& Jx, auto& Jz)
		{
			f = 0.;
			Jxdot = {{0.}};
			Jx = {{0.}};
		};

		blaze::StaticVector<Real, NX + 1, blaze::columnVector> const x0(0.);
		blaze::StaticVector<Real, NU, blaze::columnVector> const u(0.);
		blaze::StaticVector<Real, NX, blaze::columnVector> xf;

		EXPECT_THROW(irk(foo_ode, 0., 1., x0, u, xf), std::invalid_argument);
	}


	TEST(IrkTest, testInvalidDimensionU)
	{
		using Real = double;

		size_t const NX = 1, NZ = 0, NU = 1;
		ImplicitRungeKutta<Real> irk(BackwardEulerMethod(), NX, NZ, NU);

		auto foo_ode = [] (Real t, auto const& xdot, auto const& x,
			auto const& z, auto const& u, auto& f, auto& Jxdot, auto& Jx, auto& Jz)
		{
			f = 0.;
			Jxdot = {{0.}};
			Jx = {{0.}};
		};

		blaze::StaticVector<Real, NX, blaze::columnVector> const x0(0.);
		blaze::StaticVector<Real, NU + 1, blaze::columnVector> const u(0.);
		blaze::StaticVector<Real, NX, blaze::columnVector> xf;

		EXPECT_THROW(irk(foo_ode, 0., 1., x0, u, xf), std::invalid_argument);
	}
}