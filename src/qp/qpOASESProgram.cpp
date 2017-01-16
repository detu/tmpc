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

qpOASESProgram::qpOASESProgram(std::vector<QpSize> const& sz)
:	qpOASESProgram(sz, totalNumVariables(sz), totalNumConstraints(sz))
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
, 	_A(Matrix::Constant(nc, nx, Scalar{0}))
, 	_lbA(Vector::Constant(nc, sNaN))
, 	_ubA(Vector::Constant(nc, sNaN))
{
	if (sz.size() < 1)
		throw std::invalid_argument("qpOASESProgram(): at least 1 QpSize must be specified");

	nT_ = sz.size() - 1;

	// Init the -I blocks in the A matrix

	// (i, j) = top left corner of the current AB block
	size_type i = 0;
	size_type j = 0;

	for (auto sz = size_.begin(); sz + 1 < size_.end(); ++sz)
	{
		// Move (i, j) one column right from the top left corner of the current AB block,
		// which is the top left corner of the -I block.
		j += sz->nx() + sz->nu();

		// Size of the -I block
		auto const nx_next = (sz + 1)->nx();

		// Assign the -I block
		_A.block(i, j, nx_next, nx_next) = -Matrix::Identity(nx_next, nx_next);

		// Move (i, j) to the top left corner of the next AB block.
		i += nx_next + sz->nc();
	}

	// Init stage objects.
	InitStages();
}

qpOASESProgram::qpOASESProgram(qpOASESProgram const& rhs)
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

void qpOASESProgram::InitStages()
{
	assert(stage_.empty());
	stage_.reserve(size_.size());

	std::size_t const stride = nx();
	Scalar * pH = _H.data();
	Scalar * pg = _g.data();
	Scalar * plb = _lb.data();
	Scalar * pub = _ub.data();
	Scalar * pA = _A.data();
	Scalar * plbA = _lbA.data();
	Scalar * pubA = _ubA.data();

	for (auto sz = size_.cbegin(); sz != size_.cend(); ++sz)
	{
		auto const sz_next = sz + 1;
		auto const nx_next = sz_next != size_.end() ? sz_next->nx() : 0;
		stage_.emplace_back(*sz, nx_next, stride, pH, pg, plb, pub, pA, plbA, pubA);

		pH += (sz->nx() + sz->nu()) * stride + sz->nx() + sz->nu();
		pg += sz->nx() + sz->nu();
		plb += sz->nx() + sz->nu();
		pub += sz->nx() + sz->nu();
		pA += (sz->nx() + sz->nc()) * stride;
		plbA += (sz + 1)->nx() + sz->nc();
		pubA += (sz + 1)->nx() + sz->nc();
	}
}

qpOASESProgram::Stage::Stage(QpSize const& sz, std::size_t nx_next, std::size_t stride,
		Scalar * H, Scalar * g,	Scalar * lb, Scalar * ub, Scalar * A, Scalar * lbA, Scalar * ubA)
:	size_(sz)
,	Q_(H, sz.nx(), sz.nx(), Eigen::OuterStride<>(stride))
,	R_(H + sz.nx() * stride + sz.nx(), sz.nu(), sz.nu(), Eigen::OuterStride<>(stride))
,	S_(H + sz.nx(), sz.nx(), sz.nu(), Eigen::OuterStride<>(stride))
,	ST_(H + sz.nx() * stride, sz.nu(), sz.nx(), Eigen::OuterStride<>(stride))
,	q_(g, sz.nx())
,	r_(g + sz.nx(), sz.nu())
,	lbx_(lb, sz.nx())
,	ubx_(ub, sz.nx())
,	lbu_(lb + sz.nx(), sz.nu())
,	ubu_(ub + sz.nx(), sz.nu())
,	A_(A, nx_next, sz.nx(), Eigen::OuterStride<>(stride))
,	B_(A + sz.nx(), nx_next, sz.nu(), Eigen::OuterStride<>(stride))
,	lbb_(lbA, nx_next)
,	ubb_(ubA, nx_next)
,	C_(A + stride * nx_next, sz.nc(), sz.nx(), Eigen::OuterStride<>(stride))
,	D_(A + stride * nx_next + sz.nx(), sz.nc(), sz.nu(), Eigen::OuterStride<>(stride))
,	lbd_(lbA + sz.nx(), sz.nc())
,	ubd_(ubA + sz.nx(), sz.nc())
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
