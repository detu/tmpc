/*
 * CyberMotionOCP.h
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#pragma once

#include <OptimalControlProblem.hpp>

#include "CyberMotion.hpp"
#include "CasADiGeneratedFunction.hpp"

#include <Eigen/Dense>

#include <vector>

namespace mpmc
{
	class CyberMotionOCP : public camels::OptimalControlProblem<CyberMotionOCP, 2 * CyberMotion::numberOfAxes, CyberMotion::numberOfAxes>
	{
	public:
		static unsigned constexpr NY = 9;
		typedef Eigen::Matrix<Scalar, NY, 1> OutputVector;
		typedef Eigen::Matrix<Scalar, NY, NW> OutputJacobianMatrix;

		CyberMotionOCP(unsigned Nt);

		void ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac) const;
		void LagrangeTerm(unsigned i, StateInputVector const& z, StateInputVector& g, LagrangeHessianMatrix& H) const;

		StateVector const& getStateMin() const;
		StateVector const& getStateMax() const;
		InputVector const& getInputMin() const;
		InputVector const& getInputMax() const;
		StateVector const& getDefaultState() const;

		const OutputVector& getErrorWeight() const;
		void setErrorWeight(const OutputVector& val);

		void setReference(unsigned i, const OutputVector& py_ref);

		const StateVector& getWashoutState() const;
		void setWashoutState(const StateVector& val);

		double getWashoutFactor() const { return _washoutFactor; }
		void setWashoutFactor(double val) { _washoutFactor = val; }

		void MayerTerm(const StateVector& x, StateVector& g, MayerHessianMatrix& H) const;

	private:
		void Output(unsigned t, StateInputVector const& z, OutputVector& y, OutputJacobianMatrix& jac) const;

		// Private data members

		mutable CasADiGeneratedFunction _ode;
		mutable CasADiGeneratedFunction _output;

		StateVector _x_min;
		StateVector _x_max;
		InputVector _u_min;
		InputVector _u_max;
		StateVector _x_default;

		// Tracking error weights for all components of output.
		OutputVector _errorWeight;

		// Reference trajectory.
		// _yRef stores _Nt vectors of size _Ny.
		// Important: _yRef is column-major.
		std::vector<OutputVector> _yRef;

		// Washout state.
		StateVector _washoutState;

		// The more the washout factor, the more penalty for the terminal state to be far from the default (washout) position.
		double _washoutFactor;
	};
} /* namespace mpmc */

