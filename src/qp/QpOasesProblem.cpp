/*
 * qpOASESProgram.cpp
 *
 *  Created on: Apr 26, 2016
 *      Author: kotlyar
 */

#include <qp/QpOasesProblem.hpp>
#include <numeric>
#include <limits>

//#include <blaze/math/DiagonalMatrix.h>

namespace tmpc {

static auto constexpr sNaN = std::numeric_limits<QpOasesProblem::Scalar>::signaling_NaN();

static std::size_t totalNumVariables(std::vector<QpSize> const& sz)
{
	return std::accumulate(sz.begin(), sz.end(), std::size_t{0},
				[] (std::size_t n, QpSize const& sz) { return n + sz.nx() + sz.nu(); });
}

static std::size_t totalNumConstraints(std::vector<QpSize> const& sz)
{
	return std::accumulate(sz.begin(), sz.end(), std::size_t{0},
				[] (std::size_t n, QpSize const& sz) { return n + sz.nc(); })
		+ std::accumulate(sz.begin() + 1, sz.end(), std::size_t{0},
				[] (std::size_t n, QpSize const& sz) { return n + sz.nx(); });
}

QpOasesProblem::QpOasesProblem(std::vector<QpSize> const& sz)
:	QpOasesProblem(sz, totalNumVariables(sz), totalNumConstraints(sz))
{

}

QpOasesProblem::QpOasesProblem(std::initializer_list<QpSize> sz)
:	QpOasesProblem(std::vector<QpSize>(sz))
{
}

QpOasesProblem::QpOasesProblem(size_type nx, size_type nc)
:	QpOasesProblem(std::vector<QpSize>(1, QpSize(nx, 0, nc)), nx, nc)
{
}

QpOasesProblem::QpOasesProblem(std::vector<QpSize> const& sz, size_type nx, size_type nc)
:	size_(sz)
,	_H(nx, nx, Scalar{0})
,	_g(nx, sNaN)
, 	_lb(nx, sNaN)
, 	_ub(nx, sNaN)
, 	_A(nc, nx, Scalar{0})
, 	_lbA(nc, sNaN)
, 	_ubA(nc, sNaN)
{
	if (sz.size() < 1)
		throw std::invalid_argument("qpOASESProgram(): at least 1 QpSize must be specified");

	nT_ = sz.size() - 1;

	// Init the -I blocks in the A matrix and 0 blocks in lbA and ubA

	// (i, j) = top left corner of the current AB block
	size_type i = 0;
	size_type j = 0;

	for (auto sz = size_.cbegin(); sz + 1 < size_.cend(); ++sz)
	{
		// Move (i, j) one column right from the top left corner of the current AB block,
		// which is the top left corner of the -I block.
		j += sz->nx() + sz->nu();

		// Size of the -I block
		auto const nx_next = (sz + 1)->nx();

		// Assign the -I block in A
		submatrix(_A, i, j, nx_next, nx_next) = -IdentityMatrix<Matrix>(nx_next);

		// Assign the 0 blocks in lbA and ubA
		subvector(_lbA, i, nx_next) = 0.;
		subvector(_ubA, i, nx_next) = 0.;

		// Move (i, j) to the top left corner of the next AB block.
		i += nx_next + sz->nc();
	}

	// Init stage objects.
	InitStages();
}

QpOasesProblem::QpOasesProblem(QpOasesProblem const& rhs)
:	size_(rhs.size_)
,	_H(rhs._H)
,	_g(rhs._g)
, 	_lb(rhs._lb)
, 	_ub(rhs._ub)
, 	_A(rhs._A)
, 	_lbA(rhs._lbA)
, 	_ubA(rhs._ubA)
{
	InitStages();
}

void QpOasesProblem::InitStages()
{
	assert(stage_.empty());
	stage_.reserve(size_.size());

	std::size_t n = 0;
	std::size_t na = 0;

	for (auto sz = size_.cbegin(); sz != size_.cend(); ++sz)
	{
		auto const nx_next = sz + 1 != size_.end() ? (sz + 1)->nx() : 0;

		auto const Q  = submatrix(_H, n           , n           , sz->nx(), sz->nx());
		auto const S  = submatrix(_H, n           , n + sz->nx(), sz->nx(), sz->nu());
		auto const ST = submatrix(_H, n + sz->nx(), n           , sz->nu(), sz->nx());
		auto const R  = submatrix(_H, n + sz->nx(), n + sz->nx(), sz->nu(), sz->nu());

		auto const q   = subvector(_g , n           , sz->nx());
		auto const r   = subvector(_g , n + sz->nx(), sz->nu());
		auto const lbx = subvector(_lb, n           , sz->nx());
		auto const ubx = subvector(_ub, n           , sz->nx());
		auto const lbu = subvector(_lb, n + sz->nx(), sz->nu());
		auto const ubu = subvector(_ub, n + sz->nx(), sz->nu());

		auto const A = submatrix(_A, na          , n           , nx_next , sz->nx());
		auto const B = submatrix(_A, na          , n + sz->nx(), nx_next , sz->nu());
		auto const C = submatrix(_A, na + nx_next, n           , sz->nc(), sz->nx());
		auto const D = submatrix(_A, na + nx_next, n + sz->nx(), sz->nc(), sz->nx());

		auto const lbb = subvector(_lbA, na, nx_next);
		auto const ubb = subvector(_ubA, na, nx_next);
		auto const lbd = subvector(_lbA, na + nx_next, sz->nc());
		auto const ubd = subvector(_ubA, na + nx_next, sz->nc());

		stage_.emplace_back(*sz, Q, R, S, ST, q, r, lbx, ubx, lbu, ubu, A, B, lbb, ubb, C, D, lbd, ubd);

		n += sz->nx() + sz->nu();
		na += nx_next + sz->nc();
	}
}

QpOasesProblem::Stage::Stage(QpSize const& sz, SubM const& Q, SubM const& R, SubM const& S, SubM const& ST, SubV const& q, SubV const& r,
		SubV const& lbx, SubV const& ubx, SubV const& lbu, SubV const& ubu, SubM const& A, SubM const& B, SubV const& lbb, SubV const& ubb,
		SubM const& C, SubM const& D, SubV const& lbd, SubV const& ubd)
:	size_(sz)
,	Q_(Q)
,	R_(R)
,	S_(S)
,	ST_(ST)
,	q_(q)
,	r_(r)
,	lbx_(lbx)
,	ubx_(ubx)
,	lbu_(lbu)
,	ubu_(ubu)
,	A_(A)
,	B_(B)
,	lbb_(lbb)
,	ubb_(ubb)
,	C_(C)
,	D_(D)
,	lbd_(lbd)
,	ubd_(ubd)
{
}

void Print_MATLAB(std::ostream& log_stream, QpOasesProblem const& qp, std::string const& var_name)
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
