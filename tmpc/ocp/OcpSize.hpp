#pragma once

#include <tmpc/SizeT.hpp>

#include <cstdlib>
#include <vector>
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
	inline std::size_t numVariables(IteratorRange const& sz)
	{
		return numVariables(sz.begin(), sz.end());
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
	 * \brief An iterator through QP stages returning stage's OcpSize.
	 * 
	 * TODO: Deprecate?
	*/
	template <typename StageIterator>
	class QpSizeIterator
	{
	public:
		// iterator traits
		using difference_type = typename StageIterator::difference_type;
		using value_type = OcpSize const;
		using pointer = OcpSize const *;
		using reference = OcpSize const &;
		using iterator_category = typename StageIterator::iterator_category;

		explicit QpSizeIterator(StageIterator const& it)
		:	it_(it)
		{
		}

		reference operator*() const
		{
			return it_->size();
		}

		pointer operator->() const
		{
			return &it_->size();
		}

		QpSizeIterator& operator++()
		{
			++it_;
			return *this;
		}

		QpSizeIterator operator++(int)
		{
			return QpSizeIterator(it_++);
		}

		reference operator[](difference_type n) const
		{
			return it_[n].size();
		}

		friend bool operator!=(QpSizeIterator const& a, QpSizeIterator const& b)
		{
			return a.it_ != b.it_;
		}

		friend bool operator==(QpSizeIterator const& a, QpSizeIterator const& b)
		{
			return a.it_ == b.it_;
		}

		friend QpSizeIterator operator+(QpSizeIterator const& a, difference_type n)
		{
			return QpSizeIterator(a.it_ + n);
		}

		friend QpSizeIterator operator+(difference_type n, QpSizeIterator const& a)
		{
			return QpSizeIterator(n + a.it_);
		}

		friend difference_type operator-(QpSizeIterator const& a, QpSizeIterator const& b)
		{
			return a.it_ - b.it_;
		}

	private:
		StageIterator it_;
	};

	/**
	* \brief Helper function to create a QpSizeIterator from an arbitrary iterator type implementing StageIterator concept.
	*/
	template <typename StageIterator>
	inline QpSizeIterator<StageIterator> qpSizeIterator(StageIterator const& it)
	{
		return QpSizeIterator<StageIterator>(it);
	}

	/**
	* \brief Helper function to create a QpSizeIterator.
	*/
	template <typename Collection>
	inline QpSizeIterator<typename Collection::const_iterator> sizeBegin(Collection const& c)
	{
		return QpSizeIterator<typename Collection::const_iterator>(c.begin());
	}

	/**
	* \brief Helper function to create a QpSizeIterator.
	*/
	template <typename Collection>
	inline QpSizeIterator<typename Collection::const_iterator> sizeEnd(Collection const& c)
	{
		return QpSizeIterator<typename Collection::const_iterator>(c.end());
	}
}
