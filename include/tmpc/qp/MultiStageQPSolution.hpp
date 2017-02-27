/*
 * MultiStageQPSolution.hpp
 *
 *      Author: kotlyar
 */

#pragma once

#include <tmpc/Matrix.hpp>
#include <tmpc/detail/NonResizableCollection.hpp>

#include <vector>
#include <limits>

namespace tmpc
{
	template <typename Scalar_>
	class SingleStageQPSolution
	{
	public:
		typedef Scalar_ Scalar;

		SingleStageQPSolution(QpSize const& sz, size_t nx1)
		// Initialize all numeric data to NaN so that if an uninitialized object
		// by mistake used in calculations is easier to detect.
		:	size_(sz)
		,	_x(sz.nx(), sNaN())
		,	_u(sz.nu(), sNaN())
		,	_lam(2 * sz.nc() + 2 * (sz.nx() + sz.nu()), sNaN())
		{
		}

		DynamicVector<Scalar> const& get_x() const
		{
			return _x;
		}

		template <typename T>
		void set_x(T const& val)
		{
			_x = val;
		}

		DynamicVector<Scalar> const& get_u() const
		{
			return _u;
		}

		template <typename T>
		void set_u(T const& val)
		{
			_u = val;
		}

		DynamicVector<Scalar> const& get_lam() const
		{
			return _lam;
		}

		QpSize const& size() const
		{
			return size_;
		}

	private:
		QpSize size_;
		DynamicVector<Scalar> _x;
		DynamicVector<Scalar> _u;
		DynamicVector<Scalar> _lam;

		static Scalar constexpr sNaN()
	 	{
			return std::numeric_limits<Scalar>::signaling_NaN();
		}
	};
	
	//
	// Provides a generic solution interface for multistage QP solvers.
	//
	template <typename Scalar_>
	class MultiStageQPSolution : public detail::NonResizableCollection<SingleStageQPSolution<Scalar_>>
	{
	public:
		typedef Scalar_ Scalar;

		template <typename InputIterator>
		MultiStageQPSolution(InputIterator sz_first, InputIterator sz_last)
		{
			auto const N = std::distance(sz_first, sz_last);
			this->elements().reserve(N);

			for (; sz_first != sz_last; ++sz_first)
				this->elements().emplace_back(*sz_first, sz_first + 1 != sz_last ? sz_first[1].nx() : 0);
		}

		MultiStageQPSolution(MultiStageQPSolution const&) = delete;
		MultiStageQPSolution(MultiStageQPSolution&& rhs) = default;
	};
}
