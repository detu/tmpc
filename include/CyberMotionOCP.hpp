/*
 * CyberMotionOCP.h
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#pragma once

#include <OptimalControlProblem.hpp>

#include "CyberMotion.hpp"
#include "casadi_interface/GeneratedFunction.hpp"

#include <Eigen/Dense>

#include <vector>

namespace mpmc
{
	class CyberMotionOCP : public camels::OptimalControlProblem<CyberMotionOCP, 2 * CyberMotion::numberOfAxes, CyberMotion::numberOfAxes>
	{
	public:
		static unsigned const NY = 9;
		static unsigned const NC  = 2 * NU; 	// Number of non-linear constraints
		static unsigned const NCT = 2 * NU; 	// Number of non-linear terminal constraints
		typedef Eigen::Matrix<Scalar, NY, 1> OutputVector;
		typedef Eigen::Matrix<Scalar, NY, NW> OutputJacobianMatrix;
		typedef Eigen::Matrix<Scalar, NC, 1> ConstraintVector;
		typedef Eigen::Matrix<Scalar, NC, NW> ConstraintJacobianMatrix;
		typedef Eigen::Matrix<Scalar, NCT, 1> TerminalConstraintVector;
		typedef Eigen::Matrix<Scalar, NCT, NX> TerminalConstraintJacobianMatrix;

		CyberMotionOCP(unsigned Nt);

		// TODO: Consider moving integration responsibility to a separate class.
		void Integrate(const StateInputVector& z, Scalar t, StateVector& x_next, ODEJacobianMatrix& J) const;
		void LagrangeTerm(unsigned i, StateInputVector const& z, StateInputVector& g, LagrangeHessianMatrix& H) const;

		StateVector const& getStateMin() const;
		StateVector const& getStateMax() const;
		StateVector const& getTerminalStateMin() const;
		StateVector const& getTerminalStateMax() const;
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

		void PathConstraints(unsigned i, const StateInputVector& z,
			ConstraintJacobianMatrix& D, ConstraintVector& d_min, ConstraintVector& d_max) const;
		void TerminalConstraints(const StateVector& x, TerminalConstraintJacobianMatrix& D,
			TerminalConstraintVector& d_min, TerminalConstraintVector& d_max) const;

	private:
		void ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac) const;
		void Output(unsigned t, StateInputVector const& z, OutputVector& y, OutputJacobianMatrix& jac) const;
		void SRConstraints(const StateVector& x, Eigen::Matrix<Scalar, NC, NX>& D,
			ConstraintVector& d_min, ConstraintVector& d_max) const;

		// Private data members

		casadi_interface::GeneratedFunction mutable _ode;
		casadi_interface::GeneratedFunction mutable _output;

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

