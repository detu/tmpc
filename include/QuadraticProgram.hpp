#pragma once

#include <Eigen/Dense>

#include <iostream>

namespace camels
{
	template<class Scalar = double, int Options = Eigen::ColMajor>
	class QuadraticProgram
	{
	public:
		typedef unsigned int size_type;
		typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Options> Matrix;
		typedef Eigen::SelfAdjointView<Matrix, Eigen::Upper> SelfAdjointView;
		typedef Eigen::SelfAdjointView<const Matrix, Eigen::Upper> ConstSelfAdjointView;
		typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
		
		QuadraticProgram(size_type nx, size_type nc)
			: _H(nx, nx), _g(nx), _A(nc, nx), _lbA(nc), _ubA(nc), _lb(nx), _ub(nx)
		{
		}

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

		Vector& lb() { return _lb; }
		const Vector& lb() const { return _lb; }

		Vector& ub() { return _ub; }
		const Vector& ub() const { return _ub; }

		void Print_MATLAB(const std::string& var_name, std::ostream& log_stream) const
		{
			using std::endl;

			log_stream << var_name << ".H = [..." << endl << H() << "];" << endl;
			log_stream << var_name << ".g = [..." << endl << g() << "];" << endl;
			log_stream << var_name << ".A = [..." << endl << A() << "];" << endl;
			log_stream << var_name << ".lbA = [..." << endl << lbA() << "];" << endl;
			log_stream << var_name << ".ubA = [..." << endl << ubA() << "];" << endl;
			log_stream << var_name << ".lb = [..." << endl << lb() << "];" << endl;
			log_stream << var_name << ".ub = [..." << endl << ub() << "];" << endl;
		}

	private:
		Matrix _H;
		Vector _g;

		Vector _lb;
		Vector _ub;

		Matrix _A;
		Vector _lbA;
		Vector _ubA;
	};
}
