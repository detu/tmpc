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

			os << var_name << ".C{" << k + 1 << "} = [..." << endl << get_C(qp, k) << "];" << endl;
			os << var_name << ".c{" << k + 1 << "} = [..." << endl << get_c(qp, k) << "];" << endl;

			os << var_name << ".D{" << k + 1 << "} = [..." << endl << qp.get_D(k) << "];" << endl;
			os << var_name << ".dMin{" << k + 1 << "} = [..." << endl << qp.get_dMin(k) << "];" << endl;
			os << var_name << ".dMax{" << k + 1 << "} = [..." << endl << qp.get_dMax(k) << "];" << endl;

			os << var_name << ".zMin{" << k + 1 << "} = [..." << endl << get_zMin(qp, k) << "];" << endl;
			os << var_name << ".zMax{" << k + 1 << "} = [..." << endl << get_zMax(qp, k) << "];" << endl;
		}

		os << var_name << ".H{" << qp.nT() + 1 << "} = [..." << endl << get_Hend(qp) << "];" << endl;
		os << var_name << ".g{" << qp.nT() + 1 << "} = [..." << endl << get_gend(qp) << "];" << endl;

		os << var_name << ".D{" << qp.nT() + 1 << "} = [..." << endl << qp.get_Dend() << "];" << endl;
		os << var_name << ".dMin{" << qp.nT() + 1 << "} = [..." << endl << qp.get_dendMin() << "];" << endl;
		os << var_name << ".dMax{" << qp.nT() + 1 << "} = [..." << endl << qp.get_dendMax() << "];" << endl;

		os << var_name << ".zMin{" << qp.nT() + 1 << "} = [..." << endl << get_zendMin(qp) << "];" << endl;
		os << var_name << ".zMax{" << qp.nT() + 1 << "} = [..." << endl << get_zendMax(qp) << "];" << endl;
	}
}
