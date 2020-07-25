#pragma once

#include <tmpc/ocp/OcpTree.hpp>

#include <ostream>

namespace tmpc
{
	/**
	 * \brief Print an OCP tree in a human-readable form.
	 */
	inline std::ostream& operator<<(std::ostream& os, OcpTree const& g)
	{
		for (auto v : vertices(g))
		{
			os << "vertex " << v << " children: [";
			for (auto u : g.children(v))
				os << u << " ";
			os << "]" << std::endl;
		}

		return os;
	}
}
