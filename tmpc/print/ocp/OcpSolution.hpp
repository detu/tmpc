#pragma once

#include <tmpc/ocp/OcpSolutionBase.hpp>

#include <ostream>

namespace tmpc
{
	/**
	 * \brief Print an OCP solution stage in a human-readable form.
	 */
	template <typename StageSolution>
	inline std::ostream& operator<<(std::ostream& os, OcpSolutionBase<StageSolution> const& solution)
	{
		os << "x = " << trans(solution.x()) << std::endl;
		os << "u = " << trans(solution.u()) << std::endl;
		os << "pi = " << trans(solution.pi()) << std::endl;
		os << "lam_lbu = " << trans(solution.lam_lbu()) << std::endl;
		os << "lam_ubu = " << trans(solution.lam_ubu()) << std::endl;
		os << "lam_lbx = " << trans(solution.lam_lbx()) << std::endl;
		os << "lam_ubx = " << trans(solution.lam_ubx()) << std::endl;
		os << "lam_lbd = " << trans(solution.lam_lbd()) << std::endl;
		os << "lam_ubd = " << trans(solution.lam_ubd()) << std::endl;

		return os;
	}
}
