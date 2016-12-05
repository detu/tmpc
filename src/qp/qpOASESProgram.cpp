/*
 * qpOASESProgram.cpp
 *
 *  Created on: Apr 26, 2016
 *      Author: kotlyar
 */

#include <tmpc/qp/qpOASESProgram.hpp>

#include <numeric>
#include <limits>

namespace tmpc {

static auto constexpr sNaN = std::numeric_limits<qpOASESProgram::Scalar>::signaling_NaN();

qpOASESProgram::qpOASESProgram(std::vector<QpSize> const& sz)
:	qpOASESProgram(
		sz,
		std::accumulate(sz.begin(), sz.end(), size_type{0},
				[] (size_type n, QpSize const& sz) { return n + sz.nx() + sz.nu(); }),
		std::accumulate(sz.begin(), sz.end(), size_type{0},
				[] (size_type n, QpSize const& sz) { return n + sz.nc(); })
	)
{

}

qpOASESProgram::qpOASESProgram(std::initializer_list<QpSize> sz)
:	qpOASESProgram(std::vector<QpSize>(sz))
{
}

qpOASESProgram::qpOASESProgram(size_type nx, size_type nc)
:	qpOASESProgram(std::vector<QpSize>(1, QpSize(nx, 0, nc)), nx, nc)
{
}

qpOASESProgram::qpOASESProgram(std::vector<QpSize> const& sz, size_type nx, size_type nc)
:	size_(sz)
,	_H(Matrix::Constant(nx, nx, Scalar{0}))
,	_g(Vector::Constant(nx, sNaN))
, 	_lb(Vector::Constant(nx, sNaN))
, 	_ub(Vector::Constant(nx, sNaN))
, 	_A(Matrix::Constant(nc, nx, sNaN))
, 	_lbA(Vector::Constant(nc, sNaN))
, 	_ubA(Vector::Constant(nc, sNaN))
{

}

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

}	// namespace tmpc
