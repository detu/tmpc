#pragma once

#include <tmpc/ocp/OcpKktValue.hpp>

#include <ostream>

namespace tmpc
{
	/**
	 * \brief Print an OCP kkt in a human-readable form.
	 */
	template <OcpKktValue Kkt>
	inline std::ostream& operator<<(std::ostream& os, Kkt const& kkt)
	{
		auto const& g = kkt.graph();

		for (auto v : vertices(g))
		{
			bool const is_branch = out_degree(v, g) > 0;

			os << "gx[" << v << "] =\n" << trans(kkt.gx(v)) << std::endl;
			if (is_branch) os << "gu[" << v << "] =\n" << trans(kkt.gu(v)) << std::endl;
			// os << "lam_lu[" << v << "] =\n" << trans(kkt.lam_lu(v)) << std::endl;
			// os << "lam_uu[" << v << "] =\n" << trans(kkt.lam_uu(v)) << std::endl;
			// os << "lam_lx[" << v << "] =\n" << trans(kkt.lam_lx(v)) << std::endl;
			// os << "lam_ux[" << v << "] =\n" << trans(kkt.lam_ux(v)) << std::endl;
			// os << "lam_ld[" << v << "] =\n" << trans(kkt.lam_ld(v)) << std::endl;
			// os << "lam_ud[" << v << "] =\n" << trans(kkt.lam_ud(v)) << std::endl;
		}

		for (auto e : edges(g))
		{
			os << "c[" << e << "] = \n"
				<< trans(kkt.c(e)) << std::endl;
		}

		return os;
	}
}
