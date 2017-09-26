/*
 * Condensing.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: kotlyar
 */

#include <tmpc/qp/Condensing.hpp>

#include <stdexcept>
#include <numeric>

namespace tmpc {

OcpSize condensedQpSize(std::vector<OcpSize> const& sz)
{
	if (sz.empty())
		throw std::invalid_argument("CondensedSize(): vector of OcpSize must be not empty");

	return OcpSize(
		sz[0].nx(),
		std::accumulate(sz.begin(), sz.end(), std::size_t{0},
			[] (std::size_t n, OcpSize const& s) { return n + s.nu(); }),
		std::accumulate(sz.begin(), sz.end(), std::size_t{0},
			[] (std::size_t n, OcpSize const& s) { return n + s.nc(); })
		+ std::accumulate(sz.begin() + 1, sz.end(), std::size_t{0},
				[] (std::size_t n, OcpSize const& s) { return n + s.nx(); })
		);
}

} // namespace tmpc
