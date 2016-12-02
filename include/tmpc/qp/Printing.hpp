#pragma once

#include "qp.hpp"

#include <ostream>

namespace tmpc
{
	/**
	 * \brief Print a multistage QP as a MATLAB script.
	 *
	 * TODO: move all QPs to a qp namespace. Then move the printing function to the QP namespace.
	 * Then if a Print*** function is called for a QP class, it will be the function from the qp namespace.
	 */
	template <typename QP>
	inline void PrintMultistageQpMatlab(std::ostream& os, QP const& qp_, const std::string& var_name)
	{
		using std::endl;

		QP const& qp = static_cast<QP const&>(qp_);

		for (unsigned k = 0; k < qp.nT(); ++k)
		{
			os << var_name << ".H{" << k + 1 << "} = [..." << endl << get_H(qp, k) << "];" << endl;
			os << var_name << ".g{" << k + 1 << "} = [..." << endl << get_g(qp, k) << "];" << endl;

			os << var_name << ".C{" << k + 1 << "} = [..." << endl << get_AB(qp, k) << "];" << endl;
			os << var_name << ".c{" << k + 1 << "} = [..." << endl << qp.get_b(k) << "];" << endl;

			os << var_name << ".D{" << k + 1 << "} = [..." << endl << get_CD(qp, k) << "];" << endl;
			os << var_name << ".dMin{" << k + 1 << "} = [..." << endl << qp.get_d_min(k) << "];" << endl;
			os << var_name << ".dMax{" << k + 1 << "} = [..." << endl << qp.get_d_max(k) << "];" << endl;

			os << var_name << ".zMin{" << k + 1 << "} = [..." << endl << get_xu_min(qp, k) << "];" << endl;
			os << var_name << ".zMax{" << k + 1 << "} = [..." << endl << get_xu_max(qp, k) << "];" << endl;
		}

		os << var_name << ".H{" << qp.nT() + 1 << "} = [..." << endl << get_Q_end(qp) << "];" << endl;
		os << var_name << ".g{" << qp.nT() + 1 << "} = [..." << endl << get_q_end(qp) << "];" << endl;

		os << var_name << ".D{" << qp.nT() + 1 << "} = [..." << endl << qp.get_C_end() << "];" << endl;
		os << var_name << ".dMin{" << qp.nT() + 1 << "} = [..." << endl << qp.get_d_end_min() << "];" << endl;
		os << var_name << ".dMax{" << qp.nT() + 1 << "} = [..." << endl << qp.get_d_end_max() << "];" << endl;

		os << var_name << ".zMin{" << qp.nT() + 1 << "} = [..." << endl << get_x_end_min(qp) << "];" << endl;
		os << var_name << ".zMax{" << qp.nT() + 1 << "} = [..." << endl << get_x_end_max(qp) << "];" << endl;
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