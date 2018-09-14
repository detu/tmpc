#pragma once

#include <tmpc/ocp/OcpSize.hpp>

#include <vector>

namespace tmpc
{
    /**
	 * \brief Vector of OcpSize corresponding to an MPC problem with given sizes.
	 * 
	 * TODO: Obsolete. Use ocpSizeNominalMpc().
	 */
	std::vector<OcpSize> mpcOcpSize(std::size_t nt, std::size_t nx, std::size_t nu, std::size_t nc, std::size_t nct);

	/**
	 * \brief Vector of OcpSize corresponding to an MPC problem with given first stage, middle stage and last stage sizes.
	 * 
	 * TODO: Obsolete. Use ocpSizeNominalMpc().
	 */
	std::vector<OcpSize> mpcOcpSize(std::size_t nt, OcpSize const& first, OcpSize const& middle, OcpSize const& last);
}