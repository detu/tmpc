/*
 * QpSize.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: kotlyar
 */

#include <tmpc/qp/QpSize.hpp>

#include <numeric>

namespace tmpc {

std::size_t numVariables(std::vector<QpSize> const& sz)
{
	return std::accumulate(sz.begin(), sz.end(), std::size_t{0},
		[] (std::size_t n, QpSize const& s) { return n + s.nx() + s.nu(); });
}

std::size_t numEqualities(std::vector<QpSize> const& sz)
{
	return sz.empty() ? 0 : std::accumulate(sz.begin() + 1, sz.end(), std::size_t{0},
		[] (std::size_t n, QpSize const& s) { return n + s.nx(); });
}

std::size_t numInequalities(std::vector<QpSize> const& sz)
{
	return std::accumulate(sz.begin(), sz.end(), std::size_t{0},
		[] (std::size_t n, QpSize const& s) { return n + s.nx() + s.nu() + s.nc(); });
}

}	// namespace tmpc
