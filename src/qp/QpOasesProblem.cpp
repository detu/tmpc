/*
 * qpOASESProgram.cpp
 *
 *  Created on: Apr 26, 2016
 *      Author: kotlyar
 */

#include <qp/QpOasesProblem.hpp>
#include <numeric>
#include <limits>

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
		_A.block(i, j, nx_next, nx_next) = -Matrix::Identity(nx_next, nx_next);

		// Assign the 0 blocks in lbA and ubA
		_lbA.segment(i, nx_next).setZero();
		_ubA.segment(i, nx_next).setZero();

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
		auto const nx_next = sz + 1 != size_.end() ? (sz + 1)->nx() : 0;
		stage_.emplace_back(*sz, nx_next, stride, pH, pg, plb, pub, pA, plbA, pubA);

		pH += (sz->nx() + sz->nu()) * stride + sz->nx() + sz->nu();
		pg += sz->nx() + sz->nu();
		plb += sz->nx() + sz->nu();
		pub += sz->nx() + sz->nu();
		pA += (nx_next + sz->nc()) * stride + sz->nx() + sz->nu();
		plbA += nx_next + sz->nc();
		pubA += nx_next + sz->nc();
	}
}

QpOasesProblem::Stage::Stage(QpSize const& sz, std::size_t nx_next, std::size_t stride,
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
,	lbd_(lbA + nx_next, sz.nc())
,	ubd_(ubA + nx_next, sz.nc())
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
