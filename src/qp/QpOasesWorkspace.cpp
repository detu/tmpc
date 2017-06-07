#include <tmpc/qp/QpOasesWorkspace.hpp>

#include <sstream>
#include <stdexcept>

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

QpOasesSolveException::QpOasesSolveException(qpOASES::returnValue code, QpOasesProblem const& qp)
:	UnsolvedQpException("qpOASES", qp),
	_code(code),
	qpOasesProblem_(qp),
	msg_(std::string(UnsolvedQpException::what()) + "\nqpOASES return code " + std::to_string(code))
{
}

void QpOasesWorkspace::solve()
{
	/* Solve the QP. */
	int nWSR = static_cast<int>(_maxWorkingSetRecalculations);
	const auto res = _hotStart ?
		_problem.hotstart(H_.data(), g_.data(), A_.data(),
				lb_.data(), ub_.data(), lbA_.data(), ubA_.data(), nWSR) :
		_problem.init    (H_.data(), g_.data(), A_.data(),
				lb_.data(), ub_.data(), lbA_.data(), ubA_.data(), nWSR);

	if (res != qpOASES::SUCCESSFUL_RETURN)
		throw QpOasesSolveException(res, qp);

	numIter_ = nWSR;
	_hotStart = true;

	/* Get solution data. */
	_problem.getPrimalSolution(primalSolution_.data());
	_problem.getDualSolution(dualSolution_.data());
}

}	// namespace tmpc
