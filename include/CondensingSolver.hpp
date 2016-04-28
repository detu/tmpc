#pragma once

#include "MultiStageQP.hpp"
#include "qpOASESProgram.hpp"

#include <qpOASES.hpp>
#include <Condensing.hpp>

#include <ostream>

namespace camels
{
	class CondensingSolver
	{
	public:
		typedef qpOASESProgram CondensedQP;

		// Manages input data of qpOASES
		typedef camels::MultiStageQP MultiStageQP;

		// Manages output data of qpOASES
		class Point;

		// Exception that can be thrown from the Solve() member function.
		class SolveException;

		typedef unsigned size_type;
		typedef Eigen::VectorXd Vector;
		typedef Eigen::VectorXd StateVector;
		typedef Eigen::VectorXd StateInputVector;

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
		
		CondensingSolver(size_type nx, size_type nu, size_type nt)
		:	CondensingSolver(MultiStageQPSize(nx, nu, 0, 0, nt))
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

	struct CondensingSolver::SolveException : public std::runtime_error
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

	class CondensingSolver::Point
	{
	public:
		typedef Eigen::Map<Eigen::VectorXd> VectorMap;
		typedef Eigen::Map<const Eigen::VectorXd> VectorConstMap;

		Point(size_type nx, size_type nu, size_type nt);

		VectorMap w(unsigned i);
		VectorConstMap w(unsigned i) const;

		void shift();
		Point& operator+=(Point const& rhs);

		size_type const nX() const noexcept { return _nx; }
		size_type const nU() const noexcept { return _nz - _nx; }
		size_type const nT() const noexcept { return _nt; }

	private:
		size_type const _nx;
		size_type const _nz;
		size_type const _nt;

		// _data stores _Nt vectors of size _Nz and 1 vector of size _Nx
		std::vector<double> _data;
	};

	inline void CondensingSolver::Solve(MultiStageQP const& msqp, Point& solution)
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


	inline CondensingSolver::Point::Point(size_type nx, size_type nu, size_type nt)
	:	_data((nx + nu) * nt + nx)
	,	_nx(nx)
	,	_nz(nx + nu)
	,	_nt(nt)
	{
	}

	inline CondensingSolver::Point::VectorMap CondensingSolver::Point::w(unsigned i)
	{
		if(!(i < _nt + 1))
			throw std::out_of_range("CondensingSolver::Point::w(): index is out of range");

		return VectorMap(_data.data() + i * _nz, i < _nt ? _nz : _nx);
	}

	inline CondensingSolver::Point::VectorConstMap CondensingSolver::Point::w(unsigned i) const
	{
		if (!(i < _nt + 1))
			throw std::out_of_range("CondensingSolver::Point::w(): index is out of range");

		return VectorConstMap(_data.data() + i * _nz, i < _nt ? _nz : _nx);
	}

	inline void CondensingSolver::Point::shift()
	{
		std::copy_n(_data.begin() + _nz, (_nt - 1) * _nz + _nx, _data.begin());
	}

	inline CondensingSolver::Point& CondensingSolver::Point::operator+=(Point const& rhs)
	{
		if (rhs.nT() != nT())
			throw std::invalid_argument("CondensingSolver::Point::operator+=(): arguments have different sizes!");

		std::transform(_data.cbegin(), _data.cend(), rhs._data.cbegin(), _data.begin(), std::plus<double>());

		return *this;
	}
}

inline std::ostream& operator<<(std::ostream& os, camels::CondensingSolver::Point const& point)
{
	for (camels::CondensingSolver::size_type i = 0; i <= point.nT(); ++i)
		os << point.w(i) << std::endl;

	return os;
}
