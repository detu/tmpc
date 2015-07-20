#pragma once

#include <Eigen/Dense>

namespace camels
{
	class QuadraticProgram
	{
	public:
		typedef unsigned int size_type;
		typedef Eigen::MatrixXd Matrix;
		typedef Eigen::VectorXd Vector;

		QuadraticProgram(size_type nx, size_type neq, size_type nineq)
			: _H(nx, nx), _g(nx), _Aineq(nineq, nx), _bineq(nineq), _Aeq(neq, nx), _beq(neq)
		{}

		size_type nx() const { return static_cast<size_type>(_H.rows()); }
		size_type neq() const { return static_cast<size_type>(_Aeq.rows()); }
		size_type nineq() const { return static_cast<size_type>(_Aineq.rows()); }

		Matrix& H() { return _H; }
		const Matrix& H() const { return _H; }
		
		Vector& g() { return _g; }
		const Vector& g() const { return _g; }
		
		Matrix& Aeq() { return _Aeq; }
		const Matrix& Aeq() const { return _Aeq; }

		Vector& beq() { return _beq; }
		const Vector& beq() const { return _beq; }

		Matrix& Aineq() { return _Aineq; }
		const Matrix& Aineq() const { return _Aineq; }

		Vector& bineq() { return _bineq; }
		const Vector& bineq() const { return _bineq; }

	private:
		Matrix _H;
		Vector _g;

		Matrix _Aeq;
		Vector _beq;

		Matrix _Aineq;
		Vector _bineq;
	};
}
