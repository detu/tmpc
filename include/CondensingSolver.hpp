#pragma once

#include <MultiStageQP.hpp>
#include <QuadraticProgram.hpp>

#include <qpOASES.hpp>

namespace camels
{
	class CondensingSolver
	{
	public:
		typedef unsigned size_type;
		typedef QuadraticProgram<double, Eigen::RowMajor> CondensedQP;
		typedef Eigen::VectorXd Vector;

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
		void Solve(const MultiStageQP& msqp);
		const Vector& getPrimalSolution() const;
		const Vector& getPrimalCondensedSolution() const;

		const CondensedQP& getCondensedQP() const { return _condensedQP; }
		bool getHotStart() const { return _hotStart; }

	private:
		const MultiStageQPSize _size;

		// Number of constraints per stage = nX() + nD().
		size_type nC() const;

		CondensedQP _condensedQP;
		Vector _primalSolution;
		Vector _primalCondensedSolution;

		bool _hotStart = false;
		qpOASES::SQProblem _problem;
	};
}
