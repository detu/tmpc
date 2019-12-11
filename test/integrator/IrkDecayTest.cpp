#include <tmpc/integrator/ImplicitRungeKutta.hpp>
#include <tmpc/integrator/BackwardEulerMethod.hpp>
#include <tmpc/integrator/GaussLegendreMethod.hpp>

#include "DecayTest.hpp"


namespace tmpc :: testing
{
	//************************
	//
	// Exponential decay test 
	//
	//************************

	TEST_F(DecayTest, testGaussLegendre2)
	{
		testIntegrate(
			ImplicitRungeKutta<Real> {GaussLegendreMethod {2}, NX, NZ, NU, NR}
		);
	}


	TEST_F(DecayTest, testGaussLegendre2Sensitivities)
	{
		testIntegrateWithSensitivities(
			ImplicitRungeKutta<Real> {GaussLegendreMethod {2}, NX, NZ, NU, NR}
		);
	}


	TEST_F(DecayTest, testGaussLegendre2LeastSquaresLagrangeTerm)
	{
		testIntegrateLeastSquaresLagrangeTerm(
			ImplicitRungeKutta<Real> {GaussLegendreMethod {2}, NX, NZ, NU, NR}
		);
	}
}
