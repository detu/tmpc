#include <tmpc/ocp/OcpSize.hpp>

#include <numeric>

namespace tmpc {

std::size_t numVariables(std::vector<OcpSize> const& sz)
{
	return numVariables(sz.begin(), sz.end());
}

std::size_t numEqualities(std::vector<OcpSize> const& sz)
{
	return numEqualities(sz.begin(), sz.end());
}

std::size_t numInequalities(std::vector<OcpSize> const& sz)
{
	return numInequalities(sz.begin(), sz.end());
}

}	// namespace tmpc
