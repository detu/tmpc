#pragma once

#include <tmpc/ocp/OcpSolution.hpp>

#include <ostream>

namespace tmpc
{
	/**
	 * \brief Print an OCP solution in a human-readable form.
	 */
	template <OcpSolution Solution>
	inline std::ostream& operator<<(std::ostream& os, Solution const& solution)
	{
		auto const& g = solution.graph();

		for (auto v : vertices(g))
		{
			bool const is_branch = out_degree(v, g) > 0;
			
			os << "x[" << v << "] =\n" << trans(solution.x(v)) << std::endl;
			os << "lam_lx[" << v << "] =\n" << trans(solution.lam_lx(v)) << std::endl;
			os << "lam_ux[" << v << "] =\n" << trans(solution.lam_ux(v)) << std::endl;
			if (is_branch) os << "u[" << v << "] =\n" << trans(solution.u(v)) << std::endl;
			if (is_branch) os << "lam_lu[" << v << "] =\n" << trans(solution.lam_lu(v)) << std::endl;
			if (is_branch) os << "lam_uu[" << v << "] =\n" << trans(solution.lam_uu(v)) << std::endl;
			os << "lam_ld[" << v << "] =\n" << trans(solution.lam_ld(v)) << std::endl;
			os << "lam_ud[" << v << "] =\n" << trans(solution.lam_ud(v)) << std::endl;
		}

		for (auto e : edges(g))
		{
			os << "pi[" << e << "] = \n"
				<< trans(solution.pi(e)) << std::endl;
		}

		return os;
	}
}
