#pragma once

#include "MultiStageQuadraticProblemBase.hpp"

#include <ostream>

namespace tmpc
{
	template<typename QP>
	inline void Print_MATLAB(std::ostream& os, MultiStageQuadraticProblemBase<QP> const& qp_, const std::string& var_name)
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
}
