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

/**
 * \brief Total number of rows in all inequalities like lbg[k] <= C[k]*x[k] + D[k]*u[k] <= ubg[k]
 */
template <typename InputIt>
std::size_t numInequalities(InputIt sz_begin, InputIt sz_end)
{
	return std::accumulate(sz_begin, sz_end, std::size_t{0},
		[] (std::size_t n, QpSize const& s) { return n + s.nc(); });
}

/**
 * \brief An iterator through QP stages returning stage's QpSize.
 */
template <typename StageIterator>
class QpSizeIterator
{
public:
	// iterator traits
	using difference_type = typename StageIterator::difference_type;
	using value_type = QpSize const;
	using pointer = QpSize const *;
	using reference = QpSize const &;
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
inline QpSizeIterator<typename Collection::iterator> sizeBegin(Collection const& c)
{
	return QpSizeIterator<typename Collection::iterator>(c.begin());
}

/**
 * \brief Helper function to create a QpSizeIterator.
 */
template <typename Collection>
inline QpSizeIterator<typename Collection::iterator> sizeEnd(Collection const& c)
{
	return QpSizeIterator<typename Collection::iterator>(c.end());
}

}	// namespace tmpc
