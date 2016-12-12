#pragma once

#include "qpOASESProgram.hpp"
#include "QpOasesSolution.hpp"
#include "UnsolvedQpException.hpp"

#include <qpOASES.hpp>

#include <ostream>

namespace tmpc {

class QpOasesSolveException : public UnsolvedQpException
{
public:
	QpOasesSolveException(qpOASES::returnValue code, qpOASESProgram const& qp);

	qpOASES::returnValue getCode() const	{ return _code;	}
	qpOASESProgram const& getQpOasesQp() const { return qpOasesProblem_; }
	char const * what() const noexcept override { return msg_.c_str(); }

private:
	qpOASES::returnValue const _code;
	qpOASESProgram const qpOasesProblem_;
	std::string const msg_;
};

/**
 * \brief QP solver using qpOASES
 */
class QpOasesSolver
{
public:
	// Problem type for CondensingSolver
	typedef qpOASESProgram Problem;

	// Solution data type
	typedef QpOasesSolution Solution;

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
