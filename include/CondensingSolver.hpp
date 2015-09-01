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

		CondensingSolver(size_type nx, size_type nu, size_type nt);

		size_type nIndep() const { return _Nx + _Nu * _Nt; }
		size_type nDep() const { return _Nx * _Nt; }
		size_type nVar() const { return _Nz * _Nt + _Nx; }

		void Condense(const MultiStageQP& msqp);
		void Solve(const MultiStageQP& msqp);
		const Vector& getPrimalSolution() const;
		const Vector& getPrimalCondensedSolution() const;

		const CondensedQP& getCondensedQP() const { return _condensedQP; }
		bool getHotStart() const { return _hotStart; }

	private:
		size_type _Nu;
		size_type _Nx;
		size_type _Nz;
		size_type _Nt;

		CondensedQP _condensedQP;
		Vector _primalSolution;
		Vector _primalCondensedSolution;

		bool _hotStart = false;
		qpOASES::SQProblem _problem;
	};
}
