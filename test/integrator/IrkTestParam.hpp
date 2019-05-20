#pragma once

#include <tmpc/integrator/ButcherTableau.hpp>

#include <string>
#include <iosfwd>


namespace tmpc :: testing
{
	/// @brief Parameters for IRK value-parameterized tests
	template <typename Real>
	struct IrkTestParam
	{
		IrkTestParam(std::string const& method_name, ButcherTableau<Real> t, Real abs_tol, Real rel_tol)
		:	methodName(method_name)
		,	tableau(t)
		,	absTol(abs_tol)
		,	relTol(rel_tol)
		{			
		}


		/// @brief Butcher tableau for the method	
		ButcherTableau<Real> tableau;

		/// @brief Absolute tolerance for checking the integrator output
		Real absTol;

		/// @brief Relative tolerance for checking the integrator output
		Real relTol;

		/// @brief Friendly name of the method to be displayed in test output
		std::string methodName;
	};


	template <typename Real>
    inline std::ostream& operator<<(std::ostream& os, IrkTestParam<Real> const& p)
	{
		return os << p.methodName << ", absTol=" << p.absTol << ", relTol=" << p.relTol;
	}
}