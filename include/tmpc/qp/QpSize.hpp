#pragma once

#include <cstdlib>

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

}
