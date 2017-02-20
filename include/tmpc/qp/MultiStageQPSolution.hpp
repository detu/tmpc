/*
 * MultiStageQPSolution.hpp
 *
 *      Author: kotlyar
 */

#pragma once

#include <tmpc/Matrix.hpp>

#include <vector>
#include <limits>

namespace tmpc
{
	//
	// Provides a generic solution interface for multistage QP solvers.
	//
	template <typename Scalar>
	class MultiStageQPSolution
	{
	public:
		typedef Scalar_ Scalar;

		class Stage
		{
		public:
			Stage(QpSize const& sz, size_t nx1)
			// Initialize all numeric data to NaN so that if an uninitialized object
			// by mistake used in calculations is easier to detect.
			:	_x(sz.nx(), sNaN)
			,	_u(sz.nu(), sNaN)
			,	_lam(2 * sz.nc() + 2 * (sz.nx() + sz.nu()), sNaN)
			{
			}

			DynamicVector<Scalar> const& get_x() const
			{
				return _x;
			}

			DynamicVector<Scalar> const& get_u() const
			{
				return _u;
			}

			DynamicVector<Scalar> const& get_lam() const
			{
				return _lam;
			}

		private:
			DynamicVector<Scalar> _x;
			DynamicVector<Scalar> _u;
			DynamicVector<Scalar> _lam;
		};

		template <typename InputIterator>
		MultiStageQPSolution(InputIterator sz_first, InputIterator sz_last)
		{
			auto const N = std::distance(sz_first, sz_last);

			_stage.reserve(N);
			_x.reserve(N);
			_u.reserve(N);
			_lam.reserve(N);

			for (; sz_first != sz_last; ++sz_first)
			{
				_stage.emplace_back(*sz_first, sz_first + 1 != sz_last ? sz_first[1].nx() : 0);
				Stage& st = _stage.back();

				_x.push_back(st.x_data());
				_u.push_back(st.u_data());
				_lam.push_back(st.lam_data());
			}
		}

		MultiStageQPSolution(MultiStageQPSolution const&) = delete;
		MultiStageQPSolution(MultiStageQPSolution&& rhs) = default;

	private:
		static Scalar constexpr sNaN = std::numeric_limits<Scalar>::signaling_NaN();

		std::vector<Stage> _stage;
	};
}
