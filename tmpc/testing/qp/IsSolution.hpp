#pragma once

#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/qp/KktResidual.hpp>

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	/// @brief Check that the solution satisfies the given QP with the specified tolerance.
	template <typename QpSol, typename Qp, typename Real>
	inline AssertionResult isSolution(QpSol const& sol, Qp const& qp, Real abs_tol)
	{
		AssertionResult res = AssertionSuccess();
		auto const& g = qp.graph();
		auto size_map = qp.size();

		OcpKktResidual<Real> kkt_res {g, size_map};
		kktResidual(sol, qp, kkt_res);

		for (auto v : graph::vertices(g))
		{
			res = approxEqual(get(kkt_res.gx(), v), blaze::ZeroVector<Real>(get(size_map, v).nx()), abs_tol, 0.);
			if (!res)
				return res;

			res = approxEqual(get(kkt_res.gu(), v), blaze::ZeroVector<Real>(get(size_map, v).nu()), abs_tol, 0.);
			if (!res)
				return res;
		}


		for (auto e : graph::edges(g))
		{
			res = approxEqual(get(kkt_res.c(), e), blaze::ZeroVector<Real>(get(size_map, target(e, g)).nx()), abs_tol, 0.);
			if (!res)
				return res;
		}

		return res;
	}
}
