#pragma once

#include <tmpc/ocp/OcpSolution.hpp>
#include <tmpc/ocp/DynamicOcpKktValue.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/qp/KktValue.hpp>

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	/// @brief Check that the solution satisfies the given QP with the specified tolerance.
	template <OcpSolution Solution, OcpQp Qp, typename Real>
	inline AssertionResult isSolution(Solution const& sol, Qp const& qp, Real abs_tol)
	{
		AssertionResult res = AssertionSuccess();
		auto const& g = qp.graph();
		auto const& size = qp.size();

		DynamicOcpKktValue<Real> kkt_res {size};
		kktValue(qp, sol, kkt_res);

		for (auto v : vertices(g))
		{
			res = approxEqual(kkt_res.gx(v), blaze::ZeroVector<Real>(size.nx(v)), abs_tol, 0.);
			if (!res)
				return res << " gx at vertex " << v;
		}


		for (auto v : g.branchVertices())
		{
			res = approxEqual(kkt_res.gu(v), blaze::ZeroVector<Real>(size.nu(v)), abs_tol, 0.);
			if (!res)
				return res << " gu at vertex " << v;
		}


		for (auto e : edges(g))
		{
			res = approxEqual(kkt_res.c(e), blaze::ZeroVector<Real>(size.nx(target(e, g))), abs_tol, 0.);
			if (!res)
				return res << " c at edge " << e;
		}

		return res;
	}
}
