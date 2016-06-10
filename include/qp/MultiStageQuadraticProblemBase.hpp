#pragma once

//#include <Eigen/Dense>

namespace camels
{
	template<typename QP>
	decltype(auto) xMin(QP& problem, unsigned i)
	{
		return problem.zMin(i).template topRows<QP::NX>();
	}

	template<typename QP>
	decltype(auto) xMin(QP const& problem, unsigned i)
	{
		return problem.zMin(i).template topRows<QP::NX>();
	}

	template<typename QP>
	decltype(auto) xMax(QP& problem, unsigned i)
	{
		return problem.zMax(i).template topRows<QP::NX>();
	}

	template<typename QP>
	decltype(auto) xMax(QP const& problem, unsigned i)
	{
		return problem.zMax(i).template topRows<QP::NX>();
	}

	template<typename QP>
	decltype(auto) uMin(QP& problem, unsigned i)
	{
		return problem.zMin(i).template bottomRows<QP::NU>();
	}

	template<typename QP>
	decltype(auto) uMin(QP const& problem, unsigned i)
	{
		return problem.zMin(i).template bottomRows<QP::NU>();
	}

	template<typename QP>
	decltype(auto) uMax(QP& problem, unsigned i)
	{
		return problem.zMax(i).template bottomRows<QP::NU>();
	}

	template<typename QP>
	decltype(auto) uMax(QP const& problem, unsigned i)
	{
		return problem.zMax(i).template bottomRows<QP::NU>();
	}

	template<typename QP>
	typename QP::size_type nIndep(QP const& qp) { return qp.nX() + qp.nU() * qp.nT(); }

	template<typename QP>
	typename QP::size_type nDep(QP const& qp) { return qp.nX() * qp.nT(); }

	template<typename QP>
	typename QP::size_type nVar(QP const& qp) { return qp.nZ() * qp.nT() + qp.nX(); }

	template<typename QP>
	typename QP::size_type nConstr(QP const& qp) { return qp.nD() * qp.nT() + qp.nDT(); }
}
