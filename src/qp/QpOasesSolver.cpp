#include <qp/QpOasesSolver.hpp>

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

QpOasesSolveException::QpOasesSolveException(qpOASES::returnValue code, qpOASESProgram const& qp)
:	UnsolvedQpException("qpOASES", qp),
	_code(code),
	qpOasesProblem_(qp),
	msg_(std::string(UnsolvedQpException::what()) + "\nqpOASES return code " + std::to_string(code))
{
}

QpOasesSolver::QpOasesSolver(std::vector<QpSize> const& sz,	qpOASES::Options const& options)
:	size_(sz)
,	_problem(numVariables(sz), numEqualities(sz) + numInequalities(sz) - numVariables(sz))
{
	_problem.setOptions(options);
}

QpOasesSolver::QpOasesSolver(std::vector<QpSize> const& sz)
:	QpOasesSolver(sz, qpOASES_DefaultOptions())
{
}

void QpOasesSolver::Solve(Problem const& qp, Solution& solution)
{
	// Check argument sizes.
	if (qp.stageSize() != size_)
		throw std::invalid_argument("QpOasesSolver::Solve(): size of MultistageQP does not match solver sizes, sorry.");

	if (solution.size() != size_)
		throw std::invalid_argument("QpOasesSolver::Solve(): size of solution Point does not match solver sizes, sorry.");

	/* Solve the QP. */
	int nWSR = static_cast<int>(_maxWorkingSetRecalculations);
	const auto res = _hotStart ?
		_problem.hotstart(qp.H_data(), qp.g_data(), qp.A_data(),
				qp.lb_data(), qp.ub_data(), qp.lbA_data(), qp.ubA_data(), nWSR) :
		_problem.init    (qp.H_data(), qp.g_data(), qp.A_data(),
				qp.lb_data(), qp.ub_data(), qp.lbA_data(), qp.ubA_data(), nWSR);

	if (res != qpOASES::SUCCESSFUL_RETURN)
		throw QpOasesSolveException(res, qp);

	solution.setNumIter(nWSR);
	_hotStart = true;

	/* Get solution data. */
	_problem.getPrimalSolution(solution.primalSolutionData());
	_problem.getDualSolution(solution.dualSolutionData());
}

}	// namespace tmpc
