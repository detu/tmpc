#include <tmpc/qp/QpOasesWorkspace.hpp>

#include <sstream>
#include <stdexcept>
#include <limits>

namespace tmpc {

static qpOASES::Options qpOASES_DefaultOptions()
{
	qpOASES::Options options;
	options.setToMPC();
	options.printLevel = qpOASES::PL_LOW;

	assert(options.ensureConsistency() == qpOASES::SUCCESSFUL_RETURN);

	return options;
}

namespace detail {
qpOASES::Options qpOASES_DefaultOptions()
{
	return ::tmpc::qpOASES_DefaultOptions();
}
}	// namespace detail

QpOasesException::QpOasesException(qpOASES::returnValue code)
:	QpSolverException("qpOASES"),
	_code(code),
	msg_(std::string(QpSolverException::what()) + "\nqpOASES return code " + std::to_string(code))
{
}

QpOasesWorkspace::Stage::Stage(Workspace& ws, QpSize const& sz, size_t n, size_t na, size_t nx_next)
:	size_(sz)
,	Q_ (submatrix(ws.H, n          , n          , sz.nx(), sz.nx()))
,	R_ (submatrix(ws.H, n + sz.nx(), n + sz.nx(), sz.nu(), sz.nu()))
,	S_ (submatrix(ws.H, n          , n + sz.nx(), sz.nx(), sz.nu()))
,	ST_(submatrix(ws.H, n + sz.nx(), n          , sz.nu(), sz.nx()))
,	q_  (subvector(ws.g , n          , sz.nx()))
,	r_  (subvector(ws.g , n + sz.nx(), sz.nu()))
,	lbx_(subvector(ws.lb, n          , sz.nx()))
,	ubx_(subvector(ws.ub, n          , sz.nx()))
,	lbu_(subvector(ws.lb, n + sz.nx(), sz.nu()))
,	ubu_(subvector(ws.ub, n + sz.nx(), sz.nu()))
,	A_(submatrix(ws.A, na          , n          , nx_next, sz.nx()))
,	B_(submatrix(ws.A, na          , n + sz.nx(), nx_next, sz.nu()))
,	lbb_(subvector(ws.lbA, na, nx_next))
,	ubb_(subvector(ws.ubA, na, nx_next))
,	C_(submatrix(ws.A, na + nx_next, n          , sz.nc(), sz.nx()))
,	D_(submatrix(ws.A, na + nx_next, n + sz.nx(), sz.nc(), sz.nu()))
,	lbd_(subvector(ws.lbA, na + nx_next, sz.nc()))
,	ubd_(subvector(ws.ubA, na + nx_next, sz.nc()))
,	x_(subvector(ws.primalSolution, n          , sz.nx()))
,	u_(subvector(ws.primalSolution, n + sz.nx(), sz.nu()))
,	lamX_(subvector(ws.dualSolution, n          , sz.nx()))
,	lamU_(subvector(ws.dualSolution, n + sz.nx(), sz.nu()))
,	lam_(subvector(ws.dualSolution, ws.primalSolution.size() + na          , sz.nc()))
,	pi_ (subvector(ws.dualSolution, ws.primalSolution.size() + na + sz.nc(), nx_next))
{
}

QpOasesWorkspace::Workspace::Workspace(size_t nx, size_t nc)
:	H(nx, nx, Scalar{0})
,	g(nx, std::numeric_limits<Scalar>::signaling_NaN())
, 	lb(nx, std::numeric_limits<Scalar>::signaling_NaN())
, 	ub(nx, std::numeric_limits<Scalar>::signaling_NaN())
, 	A(nc, nx, Scalar{0})
, 	lbA(nc, std::numeric_limits<Scalar>::signaling_NaN())
, 	ubA(nc, std::numeric_limits<Scalar>::signaling_NaN())
,	primalSolution(nx)
,	dualSolution(nx + nc)
{		
}

void QpOasesWorkspace::solve()
{
	/* Solve the QP. */
	int nWSR = static_cast<int>(_maxWorkingSetRecalculations);
	const auto res = _hotStart ?
		problem_.hotstart(ws_.H.data(), ws_.g.data(), ws_.A.data(),
				ws_.lb.data(), ws_.ub.data(), ws_.lbA.data(), ws_.ubA.data(), nWSR) :
		problem_.init    (ws_.H.data(), ws_.g.data(), ws_.A.data(),
				ws_.lb.data(), ws_.ub.data(), ws_.lbA.data(), ws_.ubA.data(), nWSR);

	if (res != qpOASES::SUCCESSFUL_RETURN)
		throw QpOasesException(res);

	numIter_ = nWSR;
	_hotStart = true;

	/* Get solution data. */
	problem_.getPrimalSolution(ws_.primalSolution.data());
	problem_.getDualSolution(ws_.dualSolution.data());
}

}	// namespace tmpc
