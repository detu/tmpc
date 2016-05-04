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
	template<class Derived, unsigned NX_, unsigned NU_, class _Scalar = double>
	class OptimalControlProblem
	{
	public:
		static unsigned const NX = NX_;
		static unsigned const NU = NU_;
		static unsigned const NW = NX + NU;

		typedef _Scalar Scalar;
		typedef Eigen::Matrix<Scalar, NX, 1> StateVector;
		typedef Eigen::Matrix<Scalar, NU, 1> InputVector;
		typedef Eigen::Matrix<Scalar, NW, 1> StateInputVector;
		typedef Eigen::Matrix<Scalar, NX, NW> ODEJacobianMatrix;
		typedef Eigen::Matrix<Scalar, NW, NW> LagrangeHessianMatrix;
		typedef Eigen::Matrix<Scalar, NX, NX> MayerHessianMatrix;

		OptimalControlProblem(unsigned nt)
		:	_Nt(nt)
		{
		}

		void LagrangeTerm(unsigned t, StateInputVector const& z, StateInputVector& grad, LagrangeHessianMatrix& hess)
		{
			derived()->LagrangeTerm(t, z, grad, hess);
		}

		unsigned getNumberOfIntervals() const { return _Nt; }

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

		// Private data members
		const unsigned _Nt;
	};
}
