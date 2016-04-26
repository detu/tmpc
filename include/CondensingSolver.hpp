#pragma once

#include "MultiStageQP.hpp"
#include "qpOASESProgram.hpp"

#include <qpOASES.hpp>

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

	struct CondensingSolverSolveException : public std::runtime_error
	{
		CondensingSolverSolveException(qpOASES::returnValue code, qpOASESProgram const& cqp);
		qpOASES::returnValue getCode() const;
		qpOASESProgram const& getCondensedQP() const;

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
}

std::ostream& operator<<(std::ostream& os, camels::CondensingSolver::Point const& point);
