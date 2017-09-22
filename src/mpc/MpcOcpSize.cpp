#include <tmpc/mpc/MpcOcpSize.hpp>

namespace tmpc 
{
	std::vector<OcpSize> mpcOcpSize(std::size_t nt,
			std::size_t nx, std::size_t nu, std::size_t nc, std::size_t nct)
	{
		std::vector<OcpSize> sz;

		sz.reserve(nt + 1);
		std::fill_n(std::back_inserter(sz), nt, OcpSize(nx, nu, nc));
		sz.emplace_back(nx, OcpSize::size_type{0}, nct);

		return sz;
	}
}
