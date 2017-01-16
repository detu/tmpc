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
	return numVariables(sz.begin(), sz.end());
}

std::size_t numEqualities(std::vector<QpSize> const& sz)
{
	return numEqualities(sz.begin(), sz.end());
}

std::size_t numInequalities(std::vector<QpSize> const& sz)
{
	return numInequalities(sz.begin(), sz.end());
}

}	// namespace tmpc
