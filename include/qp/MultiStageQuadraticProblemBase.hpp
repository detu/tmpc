#pragma once

//#include <Eigen/Dense>

#include <cstddef>

namespace tmpc
{
	template <typename Derived>
	class MultiStageQuadraticProblemBase
	{
	public:
	};

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

	template<typename QP>
	void set_zMin(MultiStageQuadraticProblemBase<QP>& qp, std::size_t i, double val)
	{
		set_zMin(static_cast<QP&>(qp), i, QP::StateInputVector::Constant(val));
	}

	template<typename QP>
	void set_zMax(MultiStageQuadraticProblemBase<QP>& qp, std::size_t i, double val)
	{
		set_zMax(static_cast<QP&>(qp), i, QP::StateInputVector::Constant(val));
	}

	template<typename QP>
	void set_zendMin(MultiStageQuadraticProblemBase<QP>& qp, double val)
	{
		set_zendMin(static_cast<QP&>(qp), QP::StateVector::Constant(val));
	}

	template<typename QP>
	void set_zendMax(MultiStageQuadraticProblemBase<QP>& qp, double val)
	{
		set_zendMax(static_cast<QP&>(qp), QP::StateVector::Constant(val));
	}
}

namespace camels
{
	using namespace tmpc;
}
