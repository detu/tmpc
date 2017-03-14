#pragma once

#include "QpOasesProblem.hpp"
#include "QpOasesSolution.hpp"
#include "UnsolvedQpException.hpp"

#include <qpOASES.hpp>

#include <ostream>

namespace tmpc {

class QpOasesSolveException : public UnsolvedQpException
{
public:
	QpOasesSolveException(qpOASES::returnValue code, QpOasesProblem const& qp);

	qpOASES::returnValue getCode() const	{ return _code;	}
	QpOasesProblem const& getQpOasesQp() const { return qpOasesProblem_; }
	char const * what() const noexcept override { return msg_.c_str(); }

private:
	qpOASES::returnValue const _code;
	QpOasesProblem const qpOasesProblem_;
	std::string const msg_;
};

namespace detail {
qpOASES::Options qpOASES_DefaultOptions();
}

/**
 * \brief QP solver using qpOASES
 */
class QpOasesSolver
{
public:
	// Problem type for CondensingSolver
	typedef QpOasesProblem Problem;

	// Solution data type
	typedef QpOasesSolution Solution;

	template <typename InputIt>
	QpOasesSolver(InputIt sz_begin, InputIt sz_end, qpOASES::Options const& options = detail::qpOASES_DefaultOptions())
	:	size_(sz_begin, sz_end)
	,	_problem(numVariables(sz_begin, sz_end), numEqualities(sz_begin, sz_end) + numInequalities(sz_begin, sz_end))
	{
		_problem.setOptions(options);
	}

	QpOasesSolver(std::vector<QpSize> const& sz);
	QpOasesSolver(std::vector<QpSize> const& sz, qpOASES::Options const& options);

	/**
	 * \brief Copy constructor.
	 *
	 * Copying is not allowed.
	 */
	QpOasesSolver(QpOasesSolver const&) = delete;

	/**
	 * \brief Move constructor.
	 *
	 * Move-construction is ok.
	 */
	QpOasesSolver(QpOasesSolver&& rhs) = default;

	QpOasesSolver& operator=(QpOasesSolver const&) = delete;
	QpOasesSolver& operator=(QpOasesSolver&&) = delete;

	std::vector<QpSize> const& size() const
	{
		return size_;
	}

	// qpOASES-specific part
	//

	bool getHotStart() const noexcept { return _hotStart; }

	// Get maximum number of working set recalculations for qpOASES
	unsigned const getMaxWorkingSetRecalculations() const noexcept { return _maxWorkingSetRecalculations; }

	// Set maximum number of working set recalculations for qpOASES
	void setMaxWorkingSetRecalculations(unsigned val) noexcept { _maxWorkingSetRecalculations = val; }

	void Solve(Problem const& msqp, Solution& solution);

private:
	// QP size
	std::vector<QpSize> size_;

	bool _hotStart = false;

	// TODO: wrap _problem into a pImpl to
	// a) Reduce dependencies
	// b) Avoid deep-copy of qpOASES::SQProblem object of move-construction of CondensingSolver.
	qpOASES::SQProblem _problem;
	unsigned _maxWorkingSetRecalculations = 1000;
};

}	// namespace tmpc
