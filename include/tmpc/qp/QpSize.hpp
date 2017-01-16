#pragma once

#include <cstdlib>
#include <vector>
#include <numeric>

namespace tmpc {

/**
 * \brief Defines sizes of a 1-stage QP.
 */
class QpSize
{
public:
	typedef std::size_t size_type;

	QpSize(size_type const nx, size_type const nu, size_type const nc)
	:	nx_(nx)
	,	nu_(nu)
	,	nc_(nc)
	{
	}

	/**
	 * \brief Number of states.
	 */
	size_type nx() const
	{
		return nx_;
	}

	/**
	 * \brief Number of inputs.
	 */
	size_type nu() const
	{
		return nu_;
	}

	/**
	 * \brief Number of constraints.
	 */
	size_type nc() const
	{
		return nc_;
	}

private:
	size_type const nx_;
	size_type const nu_;
	size_type const nc_;
};

inline bool operator==(QpSize const& a, QpSize const& b)
{
	return a.nx() == b.nx() && a.nu() == b.nu() && a.nc() == b.nc();
}

inline bool operator!=(QpSize const& a, QpSize const& b)
{
	return !(a == b);
}

std::size_t numVariables(std::vector<QpSize> const& sz);
std::size_t numEqualities(std::vector<QpSize> const& sz);
std::size_t numInequalities(std::vector<QpSize> const& sz);

template <typename InputIt>
std::size_t numVariables(InputIt sz_begin, InputIt sz_end)
{
	return std::accumulate(sz_begin, sz_end, std::size_t{0},
		[] (std::size_t n, QpSize const& s) { return n + s.nx() + s.nu(); });
}

template <typename InputIt>
std::size_t numEqualities(InputIt sz_begin, InputIt sz_end)
{
	return sz_begin == sz_end ? 0 : std::accumulate(sz_begin + 1, sz_end, std::size_t{0},
		[] (std::size_t n, QpSize const& s) { return n + s.nx(); });
}

template <typename InputIt>
std::size_t numInequalities(InputIt sz_begin, InputIt sz_end)
{
	return std::accumulate(sz_begin, sz_end, std::size_t{0},
		[] (std::size_t n, QpSize const& s) { return n + s.nx() + s.nu() + s.nc(); });
}

}	// namespace tmpc
