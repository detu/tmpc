#pragma once

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
			os << var_name << ".Q{" << k << "} = " << AsMatlab(stage.Q()) << ";" << endl << endl;
			os << var_name << ".R{" << k << "} = " << AsMatlab(stage.R()) << ";" << endl << endl;
			os << var_name << ".S{" << k << "} = " << AsMatlab(stage.S()) << ";" << endl << endl;
			os << var_name << ".q{" << k << "} = " << AsMatlab(stage.q()) << ";" << endl << endl;
			os << var_name << ".r{" << k << "} = " << AsMatlab(stage.r()) << ";" << endl << endl;

			os << var_name << ".A{" << k << "} = " << AsMatlab(stage.A()) << ";" << endl << endl;
			os << var_name << ".B{" << k << "} = " << AsMatlab(stage.B()) << ";" << endl << endl;
			os << var_name << ".b{" << k << "} = " << AsMatlab(stage.b()) << ";" << endl << endl;

			os << var_name << ".C{" << k << "} = " << AsMatlab(stage.C()) << ";" << endl << endl;
			os << var_name << ".D{" << k << "} = " << AsMatlab(stage.D()) << ";" << endl << endl;
			os << var_name << ".lbd{" << k << "} = " << AsMatlab(stage.lbd()) << ";" << endl << endl;
			os << var_name << ".ubd{" << k << "} = " << AsMatlab(stage.ubd()) << ";" << endl << endl;

			os << var_name << ".lbx{" << k << "} = " << AsMatlab(stage.lbx()) << ";" << endl << endl;
			os << var_name << ".ubx{" << k << "} = " << AsMatlab(stage.ubx()) << ";" << endl << endl;
			os << var_name << ".lbu{" << k << "} = " << AsMatlab(stage.lbu()) << ";" << endl << endl;
			os << var_name << ".ubu{" << k << "} = " << AsMatlab(stage.ubu()) << ";" << endl << endl;
			
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
}
