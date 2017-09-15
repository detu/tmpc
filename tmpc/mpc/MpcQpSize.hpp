#pragma once

#include <tmpc/qp/QpSize.hpp>

#include <vector>

namespace tmpc
{
    /**
	 * \brief Vector of QpSize corresponding to an MPC problem with given sizes.
	 */
	std::vector<QpSize> mpcQpSize(std::size_t nt, std::size_t nx, std::size_t nu, std::size_t nc, std::size_t nct);
}