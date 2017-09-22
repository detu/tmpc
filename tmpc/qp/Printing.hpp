#pragma once

#include <tmpc/ocp/OcpSolutionBase.hpp>
#include <tmpc/qp/OcpQpBase.hpp>

#include <ostream>

namespace tmpc
{
	template <typename Matrix>
	class MatlabPrintWrapper
	{
	public:
		MatlabPrintWrapper(Matrix const& val) : val_(val) {}

		void print(std::ostream& os) const
		{
			if (val_.rows() > 0 && val_.cols() > 0)
				os << "[..." << std::endl << val_ << "]";
			else
				os << "zeros(" << val_.rows() << ", " << val_.cols() << ")";
		}

	private:
		Matrix const& val_;
	};

	template <typename Matrix>
	inline MatlabPrintWrapper<Matrix> AsMatlab(Matrix const& val)
	{
		return MatlabPrintWrapper<Matrix>(val);
	}

	template <typename Matrix>
	inline std::ostream& operator<<(std::ostream& os, MatlabPrintWrapper<Matrix> const& w)
	{
		w.print(os);
		return os;
	}

	/**
	 * \brief Print a multistage QP as a MATLAB script.
	 *
	 * TODO: move all QPs to a qp namespace. Then move the printing function to the QP namespace.
	 * Then if a Print*** function is called for a QP class, it will be the function from the qp namespace.
	 */
	template <typename QP>
	inline void PrintMultistageQpMatlab(std::ostream& os, QP const& qp, const std::string& var_name)
	{
		using std::endl;

		unsigned k = 1;

		for (auto const& stage : qp)
		{
			std::string const qp_i = var_name + "(" + std::to_string(k) + ")";

			os << qp_i << ".Q = " << AsMatlab(stage.Q()) << ";" << endl << endl;
			os << qp_i << ".R = " << AsMatlab(stage.R()) << ";" << endl << endl;
			os << qp_i << ".S = " << AsMatlab(stage.S()) << ";" << endl << endl;
			os << qp_i << ".q = " << AsMatlab(stage.q()) << ";" << endl << endl;
			os << qp_i << ".r = " << AsMatlab(stage.r()) << ";" << endl << endl;

			os << qp_i << ".A = " << AsMatlab(stage.A()) << ";" << endl << endl;
			os << qp_i << ".B = " << AsMatlab(stage.B()) << ";" << endl << endl;
			os << qp_i << ".b = " << AsMatlab(stage.b()) << ";" << endl << endl;

			os << qp_i << ".C = " << AsMatlab(stage.C()) << ";" << endl << endl;
			os << qp_i << ".D = " << AsMatlab(stage.D()) << ";" << endl << endl;
			os << qp_i << ".lbd = " << AsMatlab(stage.lbd()) << ";" << endl << endl;
			os << qp_i << ".ubd = " << AsMatlab(stage.ubd()) << ";" << endl << endl;

			os << qp_i << ".lbx = " << AsMatlab(stage.lbx()) << ";" << endl << endl;
			os << qp_i << ".ubx = " << AsMatlab(stage.ubx()) << ";" << endl << endl;
			os << qp_i << ".lbu = " << AsMatlab(stage.lbu()) << ";" << endl << endl;
			os << qp_i << ".ubu = " << AsMatlab(stage.ubu()) << ";" << endl << endl;
			
			++k;
		}
	}

	/**
	 * \brief Print a multistage QP solution in a human-readable form.
	 */
	template <typename QPSolution>
	inline void PrintMultistageQpSolution(std::ostream& os, QPSolution const& solution)
	{
		for (std::size_t i = 0; i <= solution.nT(); ++i)
		{
			os << "x[" << i << "] = " << transpose(solution.get_x(i));

			os << "\tlam_x_min[" << i << "] = " << transpose(solution.get_lam_x_min(i));
			os << "\tlam_x_max[" << i << "] = " << transpose(solution.get_lam_x_max(i));

			if (i < solution.nT())
			{
				os << "\tu[" << i << "] = " << transpose(solution.get_u(i));
				os << "\tlam_u_min[" << i << "] = " << transpose(solution.get_lam_u_min(i));
				os << "\tlam_u_max[" << i << "] = " << transpose(solution.get_lam_u_max(i));
			}

			if (i < solution.nT())
			{
				os << "\tlam_d_min[" << i << "] = " << transpose(solution.get_lam_d_min(i));
				os << "\tlam_d_max[" << i << "] = " << transpose(solution.get_lam_d_max(i));
				os << "\tpi[" << i << "] = " << transpose(solution.get_pi(i));
			}
			else
			{
				os << "\tlam_d_min[" << i << "] = " << transpose(solution.get_lam_d_end_min());
				os << "\tlam_d_max[" << i << "] = " << transpose(solution.get_lam_d_end_max());
			}

			os << std::endl;
		}
	}

	/**
	 * \brief Print a QP stage in a human-readable format.
	 */
	template <typename Stage>
	inline std::ostream& operator<<(std::ostream& os, OcpQpBase<Stage> const& stage)
	{
		using std::endl;

		os << "Q = " << endl << stage.Q() << endl << endl;
		os << "R = " << endl << stage.R() << endl << endl;
		os << "S = " << endl << stage.S() << endl << endl;
		os << "q = " << endl << trans(stage.q()) << endl << endl;
		os << "r = " << endl << trans(stage.r()) << endl << endl;

		os << "A = " << endl << stage.A() << endl << endl;
		os << "B = " << endl << stage.B() << endl << endl;
		os << "b = " << endl << trans(stage.b()) << endl << endl;

		os << "C = " << endl << stage.C() << endl << endl;
		os << "D = " << endl << stage.D() << endl << endl;
		os << "lbd = " << endl << trans(stage.lbd()) << endl << endl;
		os << "ubd = " << endl << trans(stage.ubd()) << endl << endl;

		os << "lbx = " << endl << trans(stage.lbx()) << endl << endl;
		os << "ubx = " << endl << trans(stage.ubx()) << endl << endl;
		os << "lbu = " << endl << trans(stage.lbu()) << endl << endl;
		os << "ubu = " << endl << trans(stage.ubu()) << endl << endl;

		return os;
	}

	/**
	 * \brief Print a QP stage solution in a human-readable form.
	 */
	template <typename StageSolution>
	inline std::ostream& operator<<(std::ostream& os, OcpSolutionBase<StageSolution> const& solution)
	{
		os << "x = " << trans(solution.x()) << std::endl;
		os << "u = " << trans(solution.u()) << std::endl;
		os << "pi = " << trans(solution.pi()) << std::endl;
		os << "lam_lbu = " << trans(solution.lam_lbu()) << std::endl;
		os << "lam_ubu = " << trans(solution.lam_ubu()) << std::endl;
		os << "lam_lbx = " << trans(solution.lam_lbx()) << std::endl;
		os << "lam_ubx = " << trans(solution.lam_ubx()) << std::endl;
		os << "lam_lbd = " << trans(solution.lam_lbd()) << std::endl;
		os << "lam_ubd = " << trans(solution.lam_ubd()) << std::endl;

		return os;
	}
}
