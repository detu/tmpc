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

		QuadraticProgram(size_type nx, size_type nc)
			: _H(nx, nx), _g(nx), _A(nc, nx), _lbA(nc), _ubA(nc)
		{}

		size_type nx() const { return static_cast<size_type>(_H.rows()); }
		size_type nc() const { return static_cast<size_type>(_A.rows()); }

		Matrix& H() { return _H; }
		const Matrix& H() const { return _H; }
		
		Vector& g() { return _g; }
		const Vector& g() const { return _g; }
		
		Matrix& A() { return _A; }
		const Matrix& A() const { return _A; }

		Vector& lbA() { return _lbA; }
		const Vector& lbA() const { return _lbA; }

		Vector& ubA() { return _ubA; }
		const Vector& ubA() const { return _ubA; }

	private:
		Matrix _H;
		Vector _g;

		Matrix _A;
		Vector _lbA;
		Vector _ubA;
	};
}
