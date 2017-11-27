#include <tmpc/mpc/MpcOcpSize.hpp>

namespace tmpc 
{
	std::vector<OcpSize> mpcOcpSize(size_t nt, size_t nx, size_t nu, size_t nc, size_t nct)
	{
		std::vector<OcpSize> sz;
		sz.reserve(nt + 1);

		std::fill_n(std::back_inserter(sz), nt, OcpSize {nx, nu, nc});
		sz.emplace_back(nx, size_t {0}, nct);

		return sz;
	}

	std::vector<OcpSize> mpcOcpSize(std::size_t nt, OcpSize const& first, OcpSize const& middle, OcpSize const& last)
	{
		std::vector<OcpSize> sz;		
		sz.reserve(nt + 1);

		if (nt > 0)
		{
			sz.push_back(first);
			std::fill_n(std::back_inserter(sz), nt - 1, middle);
		}
		
		sz.push_back(last);

		return sz;
	}
}
