#pragma once

#include <tmpc/ocp/OcpSize.hpp>

#include <ostream>


namespace tmpc
{
	/**
	 * \brief Print size of an OCP in human-readable format.
	 */
	template <OcpSize Size>
	inline std::ostream& operator<<(std::ostream& os, Size const& sz)
	{
		for (auto v : vertices(sz.graph()))
		{
			os << "vertex " << v << " nx=" << sz.nx(v);
			if (out_degree(v, sz.graph()) > 0)
				os << " nu=" << sz.nu(v);
			os << " nc=" << sz.nc(v);
			
			os << " children: [";
			for (auto u : sz.graph().children(v))
				os << u << " ";
			os << "]" << std::endl;
		}

		return os;
	}
}
