#include <tmpc/integrator/ExplicitRungeKutta.hpp>
#include <tmpc/integrator/RungeKutta4Method.hpp>

#include "DecayTest.hpp"


namespace tmpc :: testing
{
	//************************
	//
	// Exponential decay test 
	//
	//************************

	TEST_F(DecayTest, testErk4)
	{
		testIntegrate(ExplicitRungeKutta<Real> {RungeKutta4Method(), NX, NU, NR});
	}


	TEST_F(DecayTest, testErk4Sensitivities)
	{
		testIntegrateWithSensitivities(ExplicitRungeKutta<Real> {RungeKutta4Method(), NX, NU, NR});
	}


	TEST_F(DecayTest, testErk4LeastSquaresLagrangeTerm)
	{
		testIntegrateLeastSquaresLagrangeTerm(ExplicitRungeKutta<Real> {RungeKutta4Method(), NX, NU, NR});
	}
}