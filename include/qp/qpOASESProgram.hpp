#pragma once

#include <Eigen/Dense>

#include <iostream>

namespace tmpc
{
	// Manages input data for qpOASES solver.
	// Implements concept: QuadraticProgram.
	class qpOASESProgram
	{
		// Matrix storage option for Eigen -- important!
		// Must be RowMajor, because qpOASES expects input matrices in column-major format.
		static const int Options = Eigen::RowMajor;

	public:
		typedef unsigned int size_type;
		typedef double Scalar;
		typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;
		typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
		
		qpOASESProgram(size_type nx, size_type nc)
			: _H(nx, nx), _g(nx), _A(nc, nx), _lbA(nc), _ubA(nc), _lb(nx), _ub(nx)
		{
		}

		size_type nx() const { return static_cast<size_type>(_H.rows()); }
		size_type nc() const { return static_cast<size_type>(_A.rows()); }

		//
		// Matrix and vector access functions.
		//
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

		Vector& lb() { return _lb; }
		const Vector& lb() const { return _lb; }

		Vector& ub() { return _ub; }
		const Vector& ub() const { return _ub; }

		//
		// Raw data access functions for qpOASES.
		//
		double const * H_data() const noexcept { return _H.data(); }
		double const * g_data() const noexcept { return _g.data(); }
		double const * A_data() const noexcept { return _A.data(); }
		double const * lbA_data() const noexcept { return _lbA.data(); }
		double const * ubA_data() const noexcept { return _ubA.data(); }
		double const * lb_data() const noexcept { return _lb.data(); }
		double const * ub_data() const noexcept { return _ub.data(); }

	private:
		Matrix _H;
		Vector _g;

		Vector _lb;
		Vector _ub;

		Matrix _A;
		Vector _lbA;
		Vector _ubA;
	};

	void Print_MATLAB(std::ostream& log_stream, qpOASESProgram const& qp, std::string const& var_name);
}

namespace camels
{
	using namespace tmpc;
}
