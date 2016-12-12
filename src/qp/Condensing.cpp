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

QpSize CondensedQpSize(std::vector<QpSize> const& sz)
{
	if (sz.empty())
		throw std::invalid_argument("CondensedSize(): vector of QpSize must be not empty");

	return QpSize(
		sz[0].nx(),
		std::accumulate(sz.begin(), sz.end(), std::size_t{0},
			[] (std::size_t n, QpSize const& s) { return n + s.nu(); }),
		std::accumulate(sz.begin(), sz.end(), std::size_t{0},
			[] (std::size_t n, QpSize const& s) { return n + s.nc(); })
		);
}

} // namespace tmpc
