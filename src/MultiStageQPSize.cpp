#include <MultiStageQPSize.hpp>

#include <stdexcept>

namespace camels
{
	MultiStageQPSize::MultiStageQPSize(size_type nx, size_type nu, size_type nd, size_type ndt, size_type nt) :
		_Nx(nx),
		_Nu(nu),
		_Nt(nt),
		_Nd(nd),
		_NdT(ndt)
	{
	}

	camels::MultiStageQPSize::size_type MultiStageQPSize::nZ() const
	{
		return _Nx + _Nu;
	}

	camels::MultiStageQPSize::size_type MultiStageQPSize::nVar() const
	{
		return nZ() * nT() + nX();
	}

	camels::MultiStageQPSize::size_type MultiStageQPSize::nConstr() const
	{
		return nD() * nT() + nDT();
	}

	bool operator==(const MultiStageQPSize& s1, const MultiStageQPSize& s2)
	{
		return s1.nX() == s2.nX() && s1.nU() == s2.nU() && s1.nT() == s2.nT() && s1.nD() == s2.nD() && s1.nDT() == s2.nDT();
	}

	bool operator!=(const MultiStageQPSize& s1, const MultiStageQPSize& s2)
	{
		return !(s1 == s2);
	}

}
