/*
 * qpOASESProgram.cpp
 *
 *  Created on: Apr 26, 2016
 *      Author: kotlyar
 */

#include <qp/qpOASESProgram.hpp>

namespace tmpc
{
	void Print_MATLAB(std::ostream& log_stream, qpOASESProgram const& qp, std::string const& var_name)
	{
		using std::endl;

		log_stream << var_name << ".H = [..." << endl << qp.H() << "];" << endl;
		log_stream << var_name << ".g = [..." << endl << qp.g() << "];" << endl;
		log_stream << var_name << ".A = [..." << endl << qp.A() << "];" << endl;
		log_stream << var_name << ".lbA = [..." << endl << qp.lbA() << "];" << endl;
		log_stream << var_name << ".ubA = [..." << endl << qp.ubA() << "];" << endl;
		log_stream << var_name << ".lb = [..." << endl << qp.lb() << "];" << endl;
		log_stream << var_name << ".ub = [..." << endl << qp.ub() << "];" << endl;
	}
}
