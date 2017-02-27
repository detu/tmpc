#pragma once

#include "QuadraticProblem.hpp"
#include "QpOasesProblem.hpp"
#include "QpOasesSolver.hpp"

#include <qpOASES.hpp>
#include "Condensing.hpp"
#include "MultiStageQPSolution.hpp"
#include "UnsolvedQpException.hpp"

#include <ostream>

namespace tmpc {

/**
 * \brief Condensing solver using qpOASES
 *
 * \tparam <Scalar> Scalar type
 */
template <typename Scalar_>
class CondensingSolver
{
public:
	typedef Scalar_ Scalar;

	// Nested solver
	typedef QpOasesSolver Solver;

	// Problem type for nested solver
	typedef typename Solver::Problem CondensedQP;

	// Solution type for nested solver
	typedef typename Solver::Solution CondensedSolution;

	// Problem type for CondensingSolver
	typedef QuadraticProblem<Scalar> Problem;

	// Solution type for CondensingSolver
	typedef MultiStageQPSolution<Scalar> Solution;

	// Exception that can be thrown from the Solve() member function.
	class SolveException;

	typedef std::size_t size_type;
	typedef DynamicVector<Scalar> Vector;

	template <typename InputIterator>
	CondensingSolver(InputIterator sz_first, InputIterator sz_last)
	:	condensedSize_(condensedQpSize(sz_first, sz_last))
	,	_condensedQP({condensedSize_})
	,	_condensedSolution({condensedSize_})
	,	nestedSolver_({condensedSize_})
	{
	}

	/**
	 * \brief Copy constructor.
	 *
	 * Copying is not allowed.
	 */
	CondensingSolver(CondensingSolver const&) = delete;

	/**
	 * \brief Move constructor.
	 *
	 * Move-construction is ok.
	 */
	CondensingSolver(CondensingSolver&& rhs) = default;

	CondensingSolver& operator=(CondensingSolver const&) = delete;
	CondensingSolver& operator=(CondensingSolver&&) = delete;

	const CondensedSolution& getCondensedSolution() const { return _condensedSolution;	}

	const CondensedQP& getCondensedQP() const noexcept { return _condensedQP; }

private:
	// Size of the condensed problem
	QpSize condensedSize_;

	// Input data for nested solver
	CondensedQP _condensedQP;

	// Output data from nested solver
	CondensedSolution _condensedSolution;

	// Nested solver
	QpOasesSolver nestedSolver_;

public:
	void Solve(Problem const& msqp, Solution& solution)
	{
		// Make a condensed problem.
		_condensedQP.front() = condense<Scalar>(msqp.begin(), msqp.end());

		/* Solve the condensed QP. */
		nestedSolver_.Solve(_condensedQP, _condensedSolution);

		// Calculate the solution of the multi-stage QP.
		// TODO: add function recoverSolution() to Condensing.hpp
		auto const& condensed_solution = _condensedSolution.front();
		size_t iu = 0;
		for (size_t i = 0; i < msqp.size(); ++i)
		{
			size_t const nu = solution[i].size().nu();
			solution[i].set_u(subvector(condensed_solution.get_u(), iu, nu));

			if (i == 0)
				solution[i].set_x(condensed_solution.get_x());
			else
				solution[i].set_x(msqp[i - 1].get_A() * solution[i - 1].get_x() 
					+ msqp[i - 1].get_B() * solution[i - 1].get_u() + msqp[i - 1].get_b());

			iu += nu;
		}
	}
};

}	// namespace tmpc
