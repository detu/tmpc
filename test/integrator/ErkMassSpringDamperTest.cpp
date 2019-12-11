#include <tmpc/integrator/ExplicitRungeKutta.hpp>
#include <tmpc/integrator/RungeKutta4Method.hpp>

#include "MassSpringDamperTest.hpp"


namespace tmpc :: testing
{
	//************************
	//
	// Mass-spring-damper test 
	//
	//************************

	TEST_F(MassSpringDamperTest, testErk4)
	{
		testIntegrate(ExplicitRungeKutta<Real> {RungeKutta4Method(), NX, NU, NR});
	}


	TEST_F(MassSpringDamperTest, testErk4Sensitivities)
	{
		testIntegrateWithSensitivities(ExplicitRungeKutta<Real> {RungeKutta4Method(), NX, NU, NR});
	}


	TEST_F(MassSpringDamperTest, testErk4LeastSquaresLagrangeTerm)
	{
		testIntegrateLeastSquaresLagrangeTerm(ExplicitRungeKutta<Real> {RungeKutta4Method(), NX, NU, NR});
	}
}