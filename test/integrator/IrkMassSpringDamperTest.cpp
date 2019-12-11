#include <tmpc/integrator/ImplicitRungeKutta.hpp>
#include <tmpc/integrator/BackwardEulerMethod.hpp>
#include <tmpc/integrator/GaussLegendreMethod.hpp>

#include "MassSpringDamperTest.hpp"


namespace tmpc :: testing
{
	//************************
	//
	// Mass-spring-damper test 
	//
	//************************

	TEST_F(MassSpringDamperTest, testGaussLegendre2)
	{
		testIntegrate(
			ImplicitRungeKutta<Real> {GaussLegendreMethod {2}, NX, NZ, NU, NR}
		);
	}


	TEST_F(MassSpringDamperTest, testGaussLegendre2Sensitivities)
	{
		testIntegrateWithSensitivities(
			ImplicitRungeKutta<Real> {GaussLegendreMethod {2}, NX, NZ, NU, NR}
		);
	}


	TEST_F(MassSpringDamperTest, testGaussLegendre2LeastSquaresLagrangeTerm)
	{
		testIntegrateLeastSquaresLagrangeTerm(
			ImplicitRungeKutta<Real> {GaussLegendreMethod {2}, NX, NZ, NU, NR}
		);
	}
}