#pragma once

#include <ostream>
#include <string>

#include <Eigen/Dense>

namespace hpmpc_problem_export
{
	template <typename Formatter>
	void print_c_order_d_ip_ocp_hard_tv(Formatter& f,
										int k_max, double mu0, double mu_tol,
										int N, int const *nx, int const *nu, int const *nb, int const *ng,
										int warm_start,
										double const * const *A, double const * const *B, double const * const *b,
										double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r,
										double const * const *lb, double const * const *ub,
										double const * const *C, double const * const *D, double const * const *lg, double const * const *ug,
										double const * const *x, double const * const *u)
	{
		f.printStructure(N, nx, nu, nb, ng);

		f.printFunctionHeader();
		f.print("k_max", k_max);
		f.print("mu0", mu0);
		f.print("mu_tol", mu_tol);
		f.print("N", N);
		f.print("nx", nx, N + 1);
		f.print("nu", nu, N);
		f.print("nb", nb, N + 1);
		f.print("ng", ng, N + 1);
		f.print("warm_start", warm_start);

		for (int i = 0; i <= N; ++i)
		{
			if (i < N) f.print("A", i, A[i], nx[i + 1], nx[i]);
			if (i < N) f.print("B", i, B[i], nx[i + 1], nu[i]);
			if (i < N) f.print("b", i, b[i], nx[i + 1]);
			f.print("Q", i, Q[i], nx[i], nx[i]);
			if (i < N) f.print("S", i, S[i], nu[i], nx[i]);
			if (i < N) f.print("R", i, R[i], nu[i], nu[i]);
			f.print("q", i, q[i], nx[i]);
			if (i < N) f.print("r", i, r[i], nu[i]);
			f.print("lb", i, lb[i], nb[i]);
			f.print("ub", i, ub[i], nb[i]);
			f.print("C", i, C[i], ng[i], nx[i]);
			if (i < N) f.print("D", i, D[i], ng[i], nu[i]);
			f.print("lg", i, lg[i], ng[i]);
			f.print("ug", i, ug[i], ng[i]);
			f.print("x", i, x[i], nx[i]);
			if (i < N) f.print("u", i, u[i], nu[i]);
		}
		f.printFunctionFooter();
	}

	class MATLABFormatter
	{
	public:
		MATLABFormatter(std::ostream& os, std::string const& name_prefix = "")
		:	_os(os)
		,	_prefix(name_prefix)
		{}

		void printStructure(int N, int const *nx, int const *nu, int const *nb, int const *ng) const
		{
		}

		void printFunctionHeader() const {}
		void printFunctionFooter() const {}

		template <typename T>
		void print(std::string const& name, T const& val)
		{
			_os << assignment(name) << val << eol();
		}

		template <typename T>
		void print(std::string const& name, T const * val, std::size_t n)
		{
			_os << assignment(name);
			print(Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> const>(val, n));
			_os << eol();
		}

		template <typename T>
		void print(std::string const& name, std::size_t i, T const * val, std::size_t n)
		{
			_os << assignment(name, i);
			print(Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> const>(val, n));
			_os << eol();
		}

		template <typename T>
		void print(std::string const& name, std::size_t i, T const * val, std::size_t m, std::size_t n)
		{
			_os << assignment(name, i);
			print(Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const>(val, m, n));
			_os << eol();
		}

	private:
		template <typename Matrix>
		void print(Eigen::MatrixBase<Matrix> const& m)
		{
			if (m.rows() > 0 && m.cols() > 0)
			{
				_os << "[";

				if (m.rows() > 0)
					_os << "..." << std::endl;

				_os	<< m << "]";
			}
			else
			{
				_os << "zeros(" << m.rows() << ", " << m.cols() << ")";
			}
		}

		std::string assignment(std::string const& name) const
		{
			return _prefix + name + " = ";
		}

		std::string assignment(std::string const& name, std::size_t i) const
		{
			return _prefix + name + "{" + std::to_string(i + 1) + "} = ";
		}

		std::string eol() const { return ";\n"; }

		std::ostream& _os;
		std::string _prefix;
	};

	class CFormatter
	{
	public:
		CFormatter(std::ostream& os)
		:	_os(os)
		,	_varName("qp")
		{}

		void printStructure(int N, int const *nx, int const *nu, int const *nb, int const *ng) const;
		void printFunctionHeader() const;
		void printFunctionFooter() const;

		template <typename T>
		void print(std::string const& name, T const& val)
		{
			_os << field(name) << " = " << val << ";" << std::endl;
		}

		template <typename T>
		void print(std::string const& name, T const * val, std::size_t n)
		{
			for (std::size_t i = 0; i < n; ++i)
				_os << field(name, i) << " = " << val[i] << "; ";
			_os << std::endl;
		}

		template <typename T>
		void print(std::string const& name, std::size_t k, T const * val, std::size_t n)
		{
			std::string const name_k = name + std::to_string(k);

			if (n > 0)
			{
				for (std::size_t j = 0; j < n; ++j)
					_os << field(name_k, j) << " = " << val[j] << "; ";
				_os << std::endl;
			}

			_os << field(name, k) << " = " << field(name_k) << ";" << std::endl;
		}

		template <typename T>
		void print(std::string const& name, std::size_t k, T const * val, std::size_t m, std::size_t n)
		{
			std::string const name_k = name + std::to_string(k);

			if (n > 0)
			{
				for (std::size_t i = 0; i < m; ++i)
				{
					for (std::size_t j = 0; j < n; ++j)
						_os << field(name_k, n * i + j) << " = " << val[n * i + j] << "; ";
					_os << std::endl;
				}
			}

			_os << field(name, k) << " = " << field(name_k) << ";" << std::endl;
		}

	private:
		void printArrayStructure(std::string const& type, std::string const& name, int N) const;
		void printArrayStructure(std::string const& type, std::string const& name, int M, int N) const;
		void printArrayOfArraysStructure(std::string const& type, bool const_ptr, std::string const& name, int N, int const * n) const;
		void printArrayOfArraysStructure(std::string const& type, bool const_ptr, std::string const& name, int N, int const * m, int const * n) const;

		std::string field(std::string const& name) const
		{
			return _varName + "->" + name;
		}

		std::string field(std::string const& name, std::size_t i) const
		{
			return _varName + "->" + name + "[" + std::to_string(i) + "]";
		}

		std::string eol() const { return ";\n"; }

		std::ostream& _os;
		std::string _varName;
	};
}
