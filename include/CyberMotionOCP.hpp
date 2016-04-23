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

namespace mpmc
{
	class CyberMotionOCP : public camels::OptimalControlProblem<CyberMotionOCP, 2 * CyberMotion::numberOfAxes, CyberMotion::numberOfAxes>
	{
	public:
		static unsigned constexpr NY = 9;
		typedef Eigen::Matrix<Scalar, NY, 1> OutputVector;
		typedef Eigen::Matrix<Scalar, NY, NW> OutputJacobianMatrix;

		CyberMotionOCP();

		void ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac) const;
		void Output(unsigned t, StateInputVector const& z, OutputVector& y, OutputJacobianMatrix& jac) const;

		StateVector const& getStateMin() const;
		StateVector const& getStateMax() const;
		InputVector const& getInputMin() const;
		InputVector const& getInputMax() const;
		StateVector const& getDefaultState() const;

	private:
		mutable CasADiGeneratedFunction _ode;
		mutable CasADiGeneratedFunction _output;

		StateVector _x_min;
		StateVector _x_max;
		InputVector _u_min;
		InputVector _u_max;
		StateVector _x_default;
	};
} /* namespace mpmc */

