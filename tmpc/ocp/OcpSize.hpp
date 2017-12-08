#pragma once

#include <tmpc/SizeT.hpp>
#include <tmpc/core/TransformIterator.hpp>

#include <numeric>
#include <stdexcept>


namespace tmpc 
{
	/**
	 * \brief Defines sizes of an OCP stage.
	 *
	 * \param nx number of states
	 * \param nu number of inputs
	 * \param nc number of linear inequality constraints
	 * \param ns number of soft constraints
	*/
	class OcpSize
	{
	public:
		OcpSize()
		:	nx_ {0}
		,	nu_ {0}
		,	nc_ {0}
		,	ns_ {0}
		{
		}

		OcpSize(size_t const nx, size_t const nu, size_t const nc, size_t ns = 0)
		:	nx_ {nx}
		,	nu_ {nu}
		,	nc_ {nc}
		,	ns_ {ns}
		{
			if (ns > nx + nu + nc)
				throw std::invalid_argument("Number of slack variables is bigger than nx + nu + nc");
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
		* \brief Number of constraints.
		*/
		size_t nc() const
		{
			return nc_;
		}

		/**
		* \brief Number of soft constraints.
		*/
		size_t ns() const
		{
			return ns_;
		}

	private:
		size_t nx_;
		size_t nu_;
		size_t nc_;
		size_t ns_;
	};

	
	inline bool operator==(OcpSize const& a, OcpSize const& b)
	{
		return a.nx() == b.nx() 
			&& a.nu() == b.nu() 
			&& a.nc() == b.nc() 
			&& a.ns() == b.ns();
	}

	
	inline bool operator!=(OcpSize const& a, OcpSize const& b)
	{
		return !(a == b);
	}


	template <typename InputIt>
	size_t numVariables(InputIt sz_begin, InputIt sz_end)
	{
		return std::accumulate(sz_begin, sz_end, size_t {0},
			[] (size_t n, OcpSize const& s) { return n + s.nx() + s.nu() + 2 * s.ns(); });
	}


	template <typename IteratorRange>
	inline size_t numVariables(IteratorRange const& sz)
	{
		return numVariables(sz.begin(), sz.end());
	}


	///
	/// \return Sum of nx() of all OcpSize in the given range.
	///
	template <typename InputIt>
	inline size_t sumNx(InputIt sz_begin, InputIt sz_end)
	{
		return std::accumulate(sz_begin, sz_end, size_t {0},
			[] (size_t n, OcpSize const& s) { return n + s.nx(); });
	}


	///
	/// \return Sum of nu() of all OcpSize in the given range.
	///
	template <typename InputIt>
	inline size_t sumNu(InputIt sz_begin, InputIt sz_end)
	{
		return std::accumulate(sz_begin, sz_end, size_t {0},
			[] (size_t n, OcpSize const& s) { return n + s.nu(); });
	}


	template <typename InputIt>
	size_t numEqualities(InputIt sz_begin, InputIt sz_end)
	{
		return sz_begin == sz_end ? 0 : std::accumulate(sz_begin + 1, sz_end, size_t {0},
			[] (size_t n, OcpSize const& s) { return n + s.nx(); });
	}


	template <typename IteratorRange>
	inline std::size_t numEqualities(IteratorRange const& sz)
	{
		return numEqualities(sz.begin(), sz.end());
	}
	

	/**
	* \brief Total number of rows in all inequalities like lbg[k] <= C[k]*x[k] + D[k]*u[k] <= ubg[k]
	*/
	template <typename InputIt>
	size_t numInequalities(InputIt sz_begin, InputIt sz_end)
	{
		return std::accumulate(sz_begin, sz_end, size_t{0},
			[] (size_t n, OcpSize const& s) { return n + s.nc(); });
	}	
	

	template <typename IteratorRange>
	inline std::size_t numInequalities(IteratorRange const& sz)
	{
		return numInequalities(sz.begin(), sz.end());
	}


	/**
	* \brief Helper function to create a OcpSizeIterator from an arbitrary iterator type implementing StageIterator concept.
	*/
	template <typename StageIterator>
	inline auto ocpSizeIterator(StageIterator const& it)
	{
		return make_transform_iterator(it, [] (auto const& x) { return x.size(); });
	}


	/**
	* \brief Iterator marking the beginning of the size sequence.
	*/
	template <typename Collection>
	inline auto sizeBegin(Collection const& c)
	{
		return ocpSizeIterator(c.begin());
	}
	

	/**
	* \brief Iterator marking the end of the size sequence.
	*/
	template <typename Collection>
	inline auto sizeEnd(Collection const& c)
	{
		return ocpSizeIterator(c.end());
	}


	/**
	 * Resulting OCP size after full condensing.
	 */
	template <typename InIter>
	inline OcpSize condensedOcpSize(InIter sz_begin, InIter sz_end)
	{
		if (sz_begin == sz_end)
			throw std::invalid_argument("condensedOcpSize(): OcpSize range must be not empty");

		return OcpSize {
			sz_begin->nx(),
			std::accumulate(sz_begin, sz_end, size_t{0},
				[] (size_t n, OcpSize const& s) { return n + s.nu(); }),
			std::accumulate(sz_begin, sz_end, size_t{0},
				[] (size_t n, OcpSize const& s) { return n + s.nc(); })
			+ std::accumulate(sz_begin + 1, sz_end, size_t{0},
					[] (size_t n, OcpSize const& s) { return n + s.nx(); }),
			std::accumulate(sz_begin, sz_end, size_t{0},
				[] (size_t n, OcpSize const& s) { return n + s.ns(); })
		};
	}


	inline OcpSize condensedOcpSize(std::initializer_list<OcpSize> sz)
	{
		return condensedOcpSize(sz.begin(), sz.end());
	}
}
