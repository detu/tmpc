/*
 * CableRobotOCP.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#include <CableRobotOCP.hpp>

#include "cablerobot_generated.h"

#include <limits>

namespace mpmc
{
	CableRobotOCP::CableRobotOCP(unsigned Nt)
	:	camels::OptimalControlProblem<CableRobotOCP, 13, 8>(Nt)
	,	_ode(CASADI_GENERATED_FUNCTION_INTERFACE(cablerobot_ode))
	,	_output(CASADI_GENERATED_FUNCTION_INTERFACE(cablerobot_output))
	,	_yRef(Nt)
	,	_washoutFactor(0.)
	{
		_x_default << 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 2.5;

		auto const inf = std::numeric_limits<double>::infinity();

		_x_min << -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf, -2., -2., 1.;
		_x_max <<  inf,  inf,  inf,  inf,  inf,  inf,  inf,  inf,  inf,  inf,  2.,  2., 4.;
		_u_min << 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.;
		_u_max << 9000., 9000., 9000., 9000., 9000., 9000., 9000., 9000.;

		// Initialize error weights
		_errorWeight.fill(1.);

		// Init washout state.
		_washoutState = _x_default;
	}

	void CableRobotOCP::ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac) const
	{
		_ode({z.data(), nullptr}, {xdot.data(), jac.data()});
	}

	void CableRobotOCP::Output(unsigned t, StateInputVector const& z, OutputVector& y, OutputJacobianMatrix& jac) const
	{
		_output({z.data(), nullptr}, {y.data(), jac.data()});
	}

	CableRobotOCP::StateVector const& CableRobotOCP::getStateMin() const
	{
		return _x_min;
	}

	CableRobotOCP::StateVector const& CableRobotOCP::getStateMax() const
	{
		return _x_max;
	}

	CableRobotOCP::StateVector const& CableRobotOCP::getTerminalStateMin() const
	{
		return _x_min;
	}

	CableRobotOCP::StateVector const& CableRobotOCP::getTerminalStateMax() const
	{
		return _x_max;
	}

	CableRobotOCP::InputVector const& CableRobotOCP::getInputMin() const
	{
		return _u_min;
	}

	CableRobotOCP::InputVector const& CableRobotOCP::getInputMax() const
	{
		return _u_max;
	}

	CableRobotOCP::StateVector const& CableRobotOCP::getDefaultState() const
	{
		return _x_default;
	}

	void CableRobotOCP::LagrangeTerm(unsigned i, StateInputVector const& z, StateInputVector& g, LagrangeHessianMatrix& H) const
	{
		if (i >= getNumberOfIntervals())
			throw std::out_of_range("Time index out of range in CableRobotOCP::LagrangeTerm().");

		// Output vector and derivatives.
		OutputVector y;
		OutputJacobianMatrix G;	// G = [C, D]
		Output(i, z, y, G);

		const auto W = _errorWeight.cwiseAbs2().asDiagonal();

		// H = G^T W G + \mu I
		H = G.transpose() * W * G;

		/*
		// Quadratic term corresponding to washout (penalty for final state deviating from the default position).
		H.topLeftCorner(nX(), nX()) += _washoutFactor * MatrixXd::Identity(nX(), nX());
		*/

		// g = 2 * (y_bar - y_hat)^T * W * G
		g = (y - _yRef[i]).transpose() * W * G;

		/*
		// Linear term corresponding to washout (penalty for final state deviating from the default position).
		g.topRows(nX()) += _washoutFactor * (z.topRows(nX()) - getWashoutState());
		*/
	}

	void CableRobotOCP::setErrorWeight(const OutputVector& val)
	{
		_errorWeight = val;
	}

	const CableRobotOCP::OutputVector& CableRobotOCP::getErrorWeight() const
	{
		return _errorWeight;
	}

	void CableRobotOCP::setReference(unsigned i, const OutputVector& y_ref)
	{
		_yRef.at(i) = y_ref;
	}

	void CableRobotOCP::MayerTerm(const StateVector& x, StateVector& g, MayerHessianMatrix& H) const
	{
		// Quadratic term corresponding to washout (penalty for final state deviating from the default position).
		H = _washoutFactor * MayerHessianMatrix::Identity() * getNumberOfIntervals();

		// Linear term corresponding to washout (penalty for final state deviating from the default position).
		g = _washoutFactor * (x - getWashoutState()) * getNumberOfIntervals();
	}

	CableRobotOCP::StateVector const& CableRobotOCP::getWashoutState() const
	{
		return _washoutState;
	}

	void CableRobotOCP::setWashoutState(const StateVector& val)
	{
		_washoutState = val;
	}

	void CableRobotOCP::PathConstraints(unsigned i, const StateInputVector& z,
		ConstraintJacobianMatrix& D, ConstraintVector& d_min, ConstraintVector& d_max) const
	{

	}

	void CableRobotOCP::TerminalConstraints(const StateVector& x, TerminalConstraintJacobianMatrix& D,
		TerminalConstraintVector& d_min, TerminalConstraintVector& d_max) const
	{

	}
} /* namespace mpmc */
