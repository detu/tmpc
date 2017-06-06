#include <tmpc/mpc/MpcQpSize.hpp>

namespace tmpc {

std::vector<QpSize> mpcQpSize(std::size_t nt,
		std::size_t nx, std::size_t nu, std::size_t nc, std::size_t nct)
{
	std::vector<QpSize> sz;

	sz.reserve(nt + 1);
	std::fill_n(std::back_inserter(sz), nt, QpSize(nx, nu, nc));
	sz.emplace_back(nx, QpSize::size_type{0}, nct);

	return sz;
}

}	// namespace tmpc
