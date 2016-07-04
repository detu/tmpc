#pragma once

#include "qpOASESProgram.hpp"

#include <qpOASES.hpp>
#include <qp/Condensing.hpp>
#include <qp/MultiStageQuadraticProblem.hpp>
#include <qp/MultiStageQPSolution.hpp>

#include <ostream>

namespace tmpc
{
	class qpOASESOptions
	{
	public:
		static qpOASESOptions Reliable();
		static qpOASESOptions MPC();

		operator qpOASES::Options const&() const;
		qpOASES::PrintLevel getPrintLevel() const;
		qpOASESOptions& setPrintLevel(qpOASES::PrintLevel level);

	private:
		qpOASESOptions(qpOASES::Options const& options);

		qpOASES::Options _options;
	};

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

		// Problem type for CondensingSolver
		typedef MultiStageQuadraticProblem<NX, NU, NC, NCT> Problem;

		// Solution data type
		typedef MultiStageQPSolution<NX, NU, NC, NCT> Solution;

		// Exception that can be thrown from the Solve() member function.
		class SolveException;

		typedef unsigned size_type;
		typedef Eigen::VectorXd Vector;
		typedef Eigen::VectorXd StateVector;
		typedef Eigen::VectorXd StateInputVector;

		CondensingSolver(size_type nt, qpOASESOptions const& options = qpOASESOptions::MPC().setPrintLevel(qpOASES::PL_LOW))
		:	_condensedQP(nIndep(nt), nDep(nt) + nConstr(nt))
		,	_condensedSolution(nIndep(nt))
		,	_problem(nIndep(nt), nDep(nt) + nConstr(nt))
		,	_Nt(nt)
		{
			_problem.setOptions(options);
		}

		CondensingSolver(CondensingSolver const&) = delete;

		size_type nT() const { return _Nt; }
		size_type constexpr nX() { return NX; }
		size_type constexpr nZ() { return NZ; }
		size_type constexpr nU() { return NU; }
		size_type constexpr nD() { return NC; }
		size_type constexpr nDT() {	return NCT;	}
		size_type nIndep() const { return nIndep(nT()); }
		static size_type nIndep(size_type nt) { return NX + NU * nt; }
		size_type nDep() const { return nDep(nT()); }
		static size_type nDep(size_type nt) { return NX * nt; }
		size_type nVar() const { return nVar(nT()); }
		static size_type nVar(size_type nt) { return NZ * nt + NX; }
		static size_type nConstr(size_type nt) { return NC * nt + NCT; }

		void Solve(const Problem& msqp, Solution& solution);
		const Vector& getCondensedSolution() const { return _condensedSolution;	}

		const CondensedQP& getCondensedQP() const noexcept { return _condensedQP; }
		bool getHotStart() const noexcept { return _hotStart; }

		// qpOASES-specific part
		//

		// Get maximum number of working set recalculations for qpOASES
		unsigned const getMaxWorkingSetRecalculations() const noexcept { return _maxWorkingSetRecalculations; }

		// Set maximum number of working set recalculations for qpOASES
		void setMaxWorkingSetRecalculations(unsigned val) noexcept { _maxWorkingSetRecalculations = val; }

		// Get number of working set recalculations on last call to Solve().
		unsigned getWorkingSetRecalculations() const noexcept { return static_cast<unsigned>(_nWSR); }

	private:
		// Number of time steps
		size_type _Nt;

		// Number of constraints per stage = nX() + nD().
		size_type nC() const { return nX() + nD(); }

		// Input data for qpOASES
		CondensedQP _condensedQP;

		// Output data from qpOASES
		Vector _condensedSolution;

		// Number of working set recalculations on last call to Solve().
		int _nWSR = 0;

		bool _hotStart = false;
		qpOASES::SQProblem _problem;
		unsigned _maxWorkingSetRecalculations = 1000;
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
	void CondensingSolver<NX_, NU_, NC_, NCT_>::Solve(Problem const& msqp, Solution& solution)
	{
		// Check argument sizes.
		if (msqp.nT() != nT())
			throw std::invalid_argument("CondensingSolver::Solve(): size of MultistageQP does not match solver sizes, sorry.");

		if (solution.nT() != nT())
			throw std::invalid_argument("CondensingSolver::Solve(): size of solution Point does not match solver sizes, sorry.");

		// Make a condensed problem.
		Condense(msqp, _condensedQP);

		/* Solve the condensed QP. */
		int nWSR = static_cast<int>(_maxWorkingSetRecalculations);
		const auto res = _hotStart ?
			_problem.hotstart(_condensedQP.H_data(), _condensedQP.g_data(), _condensedQP.A_data(),
					_condensedQP.lb_data(), _condensedQP.ub_data(), _condensedQP.lbA_data(), _condensedQP.ubA_data(), nWSR) :
			_problem.init    (_condensedQP.H_data(), _condensedQP.g_data(), _condensedQP.A_data(),
					_condensedQP.lb_data(), _condensedQP.ub_data(), _condensedQP.lbA_data(), _condensedQP.ubA_data(), nWSR);

		if (res != qpOASES::SUCCESSFUL_RETURN)
			throw SolveException(res, _condensedQP);

		_nWSR = nWSR;
		_hotStart = true;

		/* Get solution of the condensed QP. */
		_problem.getPrimalSolution(_condensedSolution.data());
		//problem.getDualSolution(yOpt);

		// Calculate the solution of the multi-stage QP.
		solution.set_x(0, _condensedSolution.topRows<NX>());
		for (size_type i = 0; i < nT(); ++i)
		{
			solution.set_u(i, _condensedSolution.middleRows<NU>(NX + i * NU));
			solution.set_x(i + 1, msqp.get_A(i) * solution.get_x(i) + msqp.get_B(i) * solution.get_u(i) + msqp.get_b(i));
		}
	}
}
