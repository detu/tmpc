#pragma once

#include <tmpc/Matrix.hpp>

namespace tmpc 
{
	/**
	* \brief Defines sizes of an MPC problem.
	*/
	class MpcSize
	{
	public:
		MpcSize(size_t nx, size_t nu, size_t nc, size_t nct)
		:	nx_(nx)
		,	nu_(nu)
		,	nc_(nc)
		,	nct_(nc)
		{
		}

		/**
		* \brief Number of states.
		*/
		size_t nx() const
		{
			return nx_;
		}

		/**
		* \brief Number of inputs.
		*/
		size_t nu() const
		{
			return nu_;
		}

		/**
		* \brief Number of path constraints.
		*/
		size_t nc() const
		{
			return nc_;
		}

		/**
		* \brief Number of terminal constraints.
		*/
		size_t nct() const
		{
			return nct_;
		}

	private:
		size_t const nx_;
		size_t const nu_;
		size_t const nc_;
		size_t const nct_;
	};

	inline bool operator==(MpcSize const& a, MpcSize const& b)
	{
		return a.nx() == b.nx() && a.nu() == b.nu() && a.nc() == b.nc() && a.nct() == b.nct();
	}

	inline bool operator!=(MpcSize const& a, MpcSize const& b)
	{
		return !(a == b);
	}

}	// namespace tmpc
