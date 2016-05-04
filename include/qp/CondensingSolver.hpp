#pragma once

#include "qpOASESProgram.hpp"

#include <qpOASES.hpp>
#include <qp/Condensing.hpp>
#include <qp/qpDUNESProgram.hpp>
#include <qp/qpDUNESSolution.hpp>

#include "MultiStageQPSize.hpp"

#include <ostream>

namespace camels
{
	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	class CondensingSolver
	{
	public:
		static unsigned const NX = NX_;
		static unsigned const NU = NU_;
		static unsigned const NZ = NX + NU;
		static unsigned const NC = NC_;
		static unsigned const NCT = NCT_;

		typedef qpOASESProgram CondensedQP;

		// Manages input data of qpOASES
		typedef camels::qpDUNESProgram<NX, NU, NC, NCT> MultiStageQP;

		// Solution data type
		typedef qpDUNESSolution<NX, NU> Point;

		// Exception that can be thrown from the Solve() member function.
		class SolveException;

		typedef unsigned size_type;
		typedef Eigen::VectorXd Vector;
		typedef Eigen::VectorXd StateVector;
		typedef Eigen::VectorXd StateInputVector;

		CondensingSolver(size_type nt)
		:	CondensingSolver(MultiStageQPSize(NX, NU, NC, NCT, nt))
		{
		}		

		const MultiStageQPSize& size() const { return _size; }
		size_type nT() const { return _size.nT(); }
		size_type nX() const { return _size.nX(); }
		size_type nZ() const { return _size.nZ(); }
		size_type nU() const { return _size.nU(); }
		size_type nD() const { return _size.nD(); }
		size_type nDT() const {	return _size.nDT();	}
		size_type nIndep() const { return _size.nIndep(); }
		size_type nDep() const { return _size.nDep(); }
		size_type nVar() const { return _size.nVar(); }

		void Solve(const MultiStageQP& msqp, Point& solution);
		const Vector& getCondensedSolution() const { return _condensedSolution;	}

		const CondensedQP& getCondensedQP() const noexcept { return _condensedQP; }
		bool getHotStart() const noexcept { return _hotStart; }

	private:
		CondensingSolver(const MultiStageQPSize& size) :
			_condensedQP(size.nIndep(), size.nDep() + size.nConstr()),
			_size(size),
			_condensedSolution(size.nIndep()),
			_problem(size.nIndep(), size.nDep() + size.nConstr())
		{
			qpOASES::Options options;
			options.printLevel = qpOASES::PL_LOW;
			_problem.setOptions(options);
		}

		// Private data members.
		//
		const MultiStageQPSize _size;

		// Number of constraints per stage = nX() + nD().
		size_type nC() const { return nX() + nD(); }

		// Input data for qpOASES
		CondensedQP _condensedQP;

		// Output data from qpOASES
		Vector _condensedSolution;

		bool _hotStart = false;
		qpOASES::SQProblem _problem;
	};

	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	struct CondensingSolver<NX_, NU_, NC_, NCT_>::SolveException : public std::runtime_error
	{
		SolveException(qpOASES::returnValue code, qpOASESProgram const& cqp) :
			std::runtime_error("CondensingSolver::Solve() failed. qpOASES return code " + std::to_string(code)),
			_code(code), _CondensedQP(cqp)
		{
		}

		qpOASES::returnValue getCode() const	{ return _code;	}
		qpOASESProgram const& getCondensedQP() const { return _CondensedQP; }

	private:
		qpOASES::returnValue const _code;
		qpOASESProgram const _CondensedQP;
	};

	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	void CondensingSolver<NX_, NU_, NC_, NCT_>::Solve(MultiStageQP const& msqp, Point& solution)
	{
		// Check argument sizes.
		if (!(msqp.nX() == nX() && msqp.nU() == nU() && msqp.nT() == nT()))
			throw std::invalid_argument("CondensingSolver::Solve(): size of MultistageQP does not match solver sizes, sorry.");

		if (!(solution.nX() == nX() && solution.nU() == nU() && solution.nT() == nT()))
			throw std::invalid_argument("CondensingSolver::Solve(): size of solution Point does not match solver sizes, sorry.");

		// Make a condensed problem.
		Condense(msqp, _condensedQP);

		/* Solve the condensed QP. */
		int nWSR = 1000;
		const auto res = _hotStart ?
			_problem.hotstart(_condensedQP.H_data(), _condensedQP.g_data(), _condensedQP.A_data(),
					_condensedQP.lb_data(), _condensedQP.ub_data(), _condensedQP.lbA_data(), _condensedQP.ubA_data(), nWSR) :
			_problem.init    (_condensedQP.H_data(), _condensedQP.g_data(), _condensedQP.A_data(),
					_condensedQP.lb_data(), _condensedQP.ub_data(), _condensedQP.lbA_data(), _condensedQP.ubA_data(), nWSR);

		if (res != qpOASES::SUCCESSFUL_RETURN)
			throw SolveException(res, _condensedQP);

		_hotStart = true;

		/* Get solution of the condensed QP. */
		_problem.getPrimalSolution(_condensedSolution.data());
		//problem.getDualSolution(yOpt);

		// Calculate the solution of the multi-stage QP.
		solution.w(0).topRows(nX()) = _condensedSolution.topRows(nX());
		for (size_type i = 0; i < nT(); ++i)
		{
			auto z_i = solution.w(i);
			auto x_i = z_i.topRows(nX());
			auto u_i = z_i.bottomRows(nU());
			auto x_i_plus = solution.w(i + 1).topRows(nX());

			u_i = _condensedSolution.middleRows(nX() + i * nU(), nU());
			x_i_plus = msqp.C(i) * z_i + msqp.c(i);
		}
	}
}

// TODO: failed to write a templated version of this function.
// The following gives me "unable to deduce template arguments" error:
/*
template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
inline std::ostream& operator<<(std::ostream& os, typename camels::CondensingSolver<NX_, NU_, NC_, NCT_>::Point const& point)
{
	//typedef typename camels::CondensingSolver<NX_, NU_, NC_, NCT_>::size_type size_type;
	typedef unsigned size_type;
	for (size_type i = 0; i <= point.nT(); ++i)
		os << point.w(i) << std::endl;

	return os;
}
*/
