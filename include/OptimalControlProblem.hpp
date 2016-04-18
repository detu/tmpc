/*
 * OptimalControlProblem.hpp
 *
 *  Created on: Apr 12, 2016
 *      Author: kotlyar
 */

#pragma once

#include <Eigen/Dense>

namespace camels
{
	template<class Derived, unsigned NX, unsigned NU, unsigned NP, class _Scalar = double>
	class OptimalControlProblem
	{
	public:
		static unsigned const NW = NX + NU;

		typedef _Scalar Scalar;
		typedef Eigen::Matrix<Scalar, NX, 1> StateVector;
		typedef Eigen::Matrix<Scalar, NU, 1> InputVector;
		typedef Eigen::Matrix<Scalar, NW, 1> LagrangeGradientVector;
		typedef Eigen::Matrix<Scalar, NP, 1> ParamVector;
		typedef Eigen::Matrix<Scalar, NX, NX> StateSensitivityMatrix;
		typedef Eigen::Matrix<Scalar, NX, NU> InputSensitivityMatrix;
		typedef Eigen::Matrix<Scalar, NW, NW> LagrangeHessianMatrix;

		class ODEOutput;
		class LagrangeTermOutput;

		ODEOutput ODE(unsigned t, const StateVector& x, const InputVector& u, const ParamVector& p)
		{
			return derived()->ODE(t, x, u, p);
		}

		LagrangeTermOutput LagrangeTerm(unsigned t, const StateVector& x, const InputVector& u, const ParamVector& p)
		{
			return derived()->LagrangeTerm(t, x, u, p);
		}

		/*
		OptimalControlProblem();
		virtual ~OptimalControlProblem();

		virtual void LagrangeTerm(const Eigen::VectorXd& z, unsigned i, Eigen::MatrixXd& H, Eigen::VectorXd& g) const = 0;
		virtual void MayerTerm(const Eigen::VectorXd& x, Eigen::MatrixXd& H, Eigen::VectorXd& g) const = 0;
		virtual void PathConstraints(unsigned i, const Eigen::VectorXd& x, const Eigen::VectorXd& u, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const = 0;
		virtual void TerminalConstraints(const Eigen::VectorXd& x, Eigen::MatrixXd& D, Eigen::VectorXd& d_min, Eigen::VectorXd& d_max) const = 0;
		*/

	private:
		Derived& derived()
		{
			return static_cast<Derived&>(*this);
		}

		const Derived& derived() const
		{
			return static_cast<const Derived&>(*this);
		}
	};

	// Contains ODE derivatives xdot = f(x, u) and sensitivities d(xdot)/d(x)/ d(xdot)/d(u).
	template<class Derived, unsigned NX, unsigned NU, unsigned NP, class _Scalar>
	class OptimalControlProblem<Derived, NX, NU, NP, _Scalar>::ODEOutput
	{
	public:
		ODEOutput(const StateVector& xdot, const StateSensitivityMatrix& sens_x, const InputSensitivityMatrix& sens_u)
		:	_xdot(xdot), _sens_x(sens_x), _sens_u(sens_u)
		{
		}

		const StateVector& xDot() const
		{
			return _xdot;
		}

		const StateSensitivityMatrix& sensX() const
		{
			return _sens_x;
		}

		const InputSensitivityMatrix& sensU() const
		{
			return _sens_u;
		}

	private:
		const StateVector _xdot;
		const StateSensitivityMatrix _sens_x;
		const InputSensitivityMatrix _sens_u;
	};

	// Contains gradient and Hessian of the OCP Lagrange term.
	template<class Derived, unsigned NX, unsigned NU, unsigned NP, class _Scalar>
	class OptimalControlProblem<Derived, NX, NU, NP, _Scalar>::LagrangeTermOutput
	{
	public:
		LagrangeTermOutput(const LagrangeGradientVector& gradient, const LagrangeHessianMatrix& hessian)
		:	_gradient(gradient), _Hessian(hessian)
		{
		}

		const LagrangeGradientVector& gradient() const
		{
			return _gradient;
		}

		const LagrangeHessianMatrix& Hessian() const
		{
			return _Hessian;
		}

	private:
		const LagrangeGradientVector _gradient;
		const LagrangeHessianMatrix _Hessian;
	};
}
