#pragma once

namespace camels
{
	/**
	MultiStageQPSize defines sizes of a multistage QP problem.
	*/
	class MultiStageQPSize
	{
	public:
		typedef unsigned int size_type;

		MultiStageQPSize(size_type nx, size_type nu, size_type nd, size_type ndt, size_type nt);

		size_type nT() const { return _Nt; }
		size_type nX() const { return _Nx; }
		size_type nZ() const;
		size_type nU() const { return _Nu; }
		size_type nD() const { return _Nd; }
		size_type nDT() const { return _NdT; }
		size_type nIndep() const { return _Nx + _Nu * _Nt; }
		size_type nDep() const { return _Nx * _Nt; }
		size_type nVar() const;

		// Number of path constraints for all stages + number of terminal constraints.
		size_type nConstr() const;

	private:
		const size_type _Nu;
		const size_type _Nx;
		const size_type _Nt;
		const size_type _Nd;
		const size_type _NdT;
	};

	bool operator==(const MultiStageQPSize& s1, const MultiStageQPSize& s2);
	bool operator!=(const MultiStageQPSize& s1, const MultiStageQPSize& s2);
}
