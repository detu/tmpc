/*
 * CyberMotionOCP.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#include <cms/CyberMotionOCP.hpp>

#include "cybermotion_generated.h"

namespace cms
{
	CyberMotionOCP::CyberMotionOCP(unsigned Nt)
	:	camels::OptimalControlProblem<CyberMotionOCP, 2 * CyberMotion::numberOfAxes, CyberMotion::numberOfAxes>(Nt)
	,	_ode(CASADI_GENERATED_FUNCTION_INTERFACE(cybermotion_ode))
	,	_output(CASADI_GENERATED_FUNCTION_INTERFACE(cybermotion_output))
	,	_yRef(Nt + 1)
	,	_washoutFactor(0.)
	{
		// Initialize state and input limits.

		static const double limits[NU][6] = {
			{ 0.2, 9.3, -1.47, 1.47, -1.0780, 1.0780 },
			{ -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), -1.1802, 1.1802, -1.6762, 1.6762 },
			{ -2.1447, -0.8727, -0.9749, 0.9749, -1.1973, 1.1973 },
			{ -0.7540, 1.5415, -1.1802, 1.1802, -2.1893, 2.1893 },
			{ -3.0159, 3.0159, -1.2999, 1.2999, -0.5644, 0.5644 },
			{ -0.8727, 0.8727, -1.2999, 1.2999, -1.6249, 1.6249 },
			{ -3.0159, 3.0159, -2.0525, 2.0525, -1.3170, 1.3170 },
			// Full-range position for cabin axis is { 0.1938, 1.4839, -0.2450, 0.2450, -0.9800, 0.9800 }.
			// Below, the cabin axis position is limited to both rollers on the curve.
			// This is the only case implemented in current CasADi-generated code.
			{ 0.5700, 1.1603, -0.2450, 0.2450, -0.9800, 0.9800 }
		};

		for (int i = 0; i < NU; ++i)
		{
			_x_min[i] 	   = limits[i][0];	_x_max[i] 	   = limits[i][1];
			_x_min[NU + i] = limits[i][2];	_x_max[NU + i] = limits[i][3];
			_u_min[i] 	   = limits[i][4];	_u_max[i]	   = limits[i][5];
		}

		// The "agile" position and zero velocity.
		_x_default << 4.8078, 0.1218, -1.5319, 0.4760, 0.0006, 0.1396, -0.0005, 0.7991, 0., 0., 0., 0., 0., 0., 0., 0.;

		// Initialize error weights
		_errorWeight.fill(1.);

		// Init washout state.
		_washoutState = _x_default;
	}

	void CyberMotionOCP::ODE(unsigned t, StateInputVector const& z, StateVector& xdot, ODEJacobianMatrix& jac) const
	{
		_ode({z.data(), nullptr}, {xdot.data(), jac.data()});
	}

	void CyberMotionOCP::Output(unsigned t, StateInputVector const& z, OutputVector& y, OutputJacobianMatrix& jac) const
	{
		_output({z.data(), nullptr}, {y.data(), jac.data()});
	}

	CyberMotionOCP::StateVector const& CyberMotionOCP::getStateMin() const
	{
		return _x_min;
	}

	CyberMotionOCP::StateVector const& CyberMotionOCP::getStateMax() const
	{
		return _x_max;
	}

	CyberMotionOCP::StateVector const& CyberMotionOCP::getTerminalStateMin() const
	{
		return _x_min;
	}

	CyberMotionOCP::StateVector const& CyberMotionOCP::getTerminalStateMax() const
	{
		return _x_max;
	}

	CyberMotionOCP::InputVector const& CyberMotionOCP::getInputMin() const
	{
		return _u_min;
	}

	CyberMotionOCP::InputVector const& CyberMotionOCP::getInputMax() const
	{
		return _u_max;
	}

	CyberMotionOCP::StateVector const& CyberMotionOCP::getDefaultState() const
	{
		return _x_default;
	}

	void CyberMotionOCP::LagrangeTerm(unsigned i, StateInputVector const& z, StateInputVector& g, LagrangeHessianMatrix& H) const
	{
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

	void CyberMotionOCP::setErrorWeight(const OutputVector& val)
	{
		_errorWeight = val;
	}

	const CyberMotionOCP::OutputVector& CyberMotionOCP::getErrorWeight() const
	{
		return _errorWeight;
	}

	void CyberMotionOCP::setReference(unsigned i, const OutputVector& y_ref)
	{
		_yRef.at(i) = y_ref;
	}

	void CyberMotionOCP::MayerTerm(const StateVector& x, StateVector& g, MayerHessianMatrix& H) const
	{
		// Quadratic term corresponding to washout (penalty for final state deviating from the default position).
		H = _washoutFactor * MayerHessianMatrix::Identity() * getNumberOfIntervals();

		// Linear term corresponding to washout (penalty for final state deviating from the default position).
		g = _washoutFactor * (x - getWashoutState()) * getNumberOfIntervals();
	}

	CyberMotionOCP::StateVector const& CyberMotionOCP::getWashoutState() const
	{
		return _washoutState;
	}

	void CyberMotionOCP::setWashoutState(const StateVector& val)
	{
		_washoutState = val;
	}

	void CyberMotionOCP::SRConstraints(const StateVector& x, Eigen::Matrix<Scalar, NC, NX>& D, ConstraintVector& d_min, ConstraintVector& d_max) const
	{
		const auto q = x.topRows<NU>();
		const auto v = x.bottomRows<NU>();
		const auto q_min = _x_min.topRows<NU>();
		const auto q_max = _x_max.topRows<NU>();

		typedef Eigen::Matrix<Scalar, NU, NU> MatrixUU;

		/*
		Stoppability-reachability equations:
		v^2 + 2 * (q_max - q) * a_min <= 0
		v^2 + 2 * (q_min - q) * a_max <= 0

		Linearization:
		-inf <= 2 * (-a_min * dq + v * dv) <= -v^2 - 2 * (q_max - q) * a_min
		-inf <= 2 * (-a_min * dq + v * dv) <= -v^2 - 2 * (q_min - q) * a_max
		*/
		D << -2. * MatrixUU(_u_min.asDiagonal()), 2. * MatrixUU(q.asDiagonal()),
			 -2. * MatrixUU(_u_max.asDiagonal()), 2. * MatrixUU(q.asDiagonal());

		const auto inf = InputVector::Constant(std::numeric_limits<double>::infinity());
		d_min << -inf, -inf;
		d_max << -v.cwiseAbs2() - 2. * (q_max - q).cwiseProduct(getInputMin()), -v.cwiseAbs2() - 2. * (q_min - q).cwiseProduct(getInputMax());
	}

	void CyberMotionOCP::PathConstraints(unsigned i, const StateInputVector& z,
		ConstraintJacobianMatrix& D, ConstraintVector& d_min, ConstraintVector& d_max) const
	{
		Eigen::Matrix<Scalar, NC, NX> Dx;
		SRConstraints(z.topRows<NX>(), Dx, d_min, d_max);
		D << Dx, Eigen::Matrix<Scalar, NC, NU>::Zero();
	}

	void CyberMotionOCP::TerminalConstraints(const StateVector& x, TerminalConstraintJacobianMatrix& D,
		TerminalConstraintVector& d_min, TerminalConstraintVector& d_max) const
	{
		SRConstraints(x, D, d_min, d_max);
	}

	void CyberMotionOCP::Integrate(const StateInputVector& z, Scalar const t, StateVector& x_next, ODEJacobianMatrix& J) const
	{
		static auto const I = Eigen::Matrix<double, 8, 8>::Identity();
		static auto const O = Eigen::Matrix<double, 8, 8>::Zero();

		J << I, t * I, std::pow(t, 2) / 2. * I,
			 O,     I,                   t * I;

		x_next = J * z;
	}
} /* namespace mpmc */
