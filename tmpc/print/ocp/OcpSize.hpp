#pragma once

#include <tmpc/ocp/OcpSize.hpp>

#include <ostream>

namespace tmpc
{
	/**
	 * \brief Print size of an OCP QP in a human-readable format.
	 */
	inline std::ostream& operator<<(std::ostream& os, OcpSize const& sz)
	{
		return os << "OcpSize(nx=" << sz.nx() << ", nu=" << sz.nu() << ", nc=" << sz.nc() << ")";
	}
}
