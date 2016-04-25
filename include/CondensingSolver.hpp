#pragma once

#include <MultiStageQP.hpp>
#include <QuadraticProgram.hpp>

#include <qpOASES.hpp>

#include <ostream>

namespace camels
{
	class CondensingSolver
	{
	public:
		// Manages input data of qpOASES
		typedef camels::MultiStageQP MultiStageQP;

		// Manages output data of qpOASES
		class Point;

		typedef unsigned size_type;
		typedef QuadraticProgram<double, Eigen::RowMajor> CondensedQP;
		typedef Eigen::VectorXd Vector;
		typedef Eigen::VectorXd StateVector;
		typedef Eigen::VectorXd StateInputVector;

		CondensingSolver(const MultiStageQPSize& size);
		CondensingSolver(size_type nx, size_type nu, size_type nt);

		const MultiStageQPSize& size() const;
		size_type nT() const;
		size_type nX() const;
		size_type nZ() const;
		size_type nU() const;
		size_type nD() const;
		size_type nDT() const;
		size_type nIndep() const;
		size_type nDep() const;
		size_type nVar() const;

		void Condense(const MultiStageQP& msqp);
		void Solve(const MultiStageQP& msqp, Point& solution);
		const Vector& getPrimalCondensedSolution() const;

		const CondensedQP& getCondensedQP() const noexcept { return _condensedQP; }
		bool getHotStart() const noexcept { return _hotStart; }

	private:
		const MultiStageQPSize _size;

		// Number of constraints per stage = nX() + nD().
		size_type nC() const;

		CondensedQP _condensedQP;
		Vector _primalCondensedSolution;

		bool _hotStart = false;
		qpOASES::SQProblem _problem;
	};

	struct CondensingSolverSolveException : public std::runtime_error
	{
		CondensingSolverSolveException(qpOASES::returnValue code, const CondensingSolver::CondensedQP& cqp);
		const qpOASES::returnValue getCode() const;
		const CondensingSolver::CondensedQP& getCondensedQP() const;

	private:
		const qpOASES::returnValue _code;
		const CondensingSolver::CondensedQP _CondensedQP;
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
