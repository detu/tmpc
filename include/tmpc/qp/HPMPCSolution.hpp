/*
 * qpDUNESSolution.hpp
 *
 *  Created on: May 4, 2016
 *      Author: kotlyar
 */

#pragma once

#include <tmpc/Matrix.hpp>

#include <vector>
#include <limits>

namespace tmpc
{
	//
	// Provides solution interface for the HPMPC solver.
	//
	template <typename Scalar_>
	class HPMPCSolution
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
			,	_pi(nx1, sNaN)
			,	_lam(2 * sz.nc() + 2 * (sz.nx() + sz.nu()), sNaN)
			,	_t(2 * sz.nc() + 2 * (sz.nx() + sz.nu()), sNaN)
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

			DynamicVector<Scalar> const& get_pi() const
			{
				return _pi;
			}

			DynamicVector<Scalar> const& get_lam() const
			{
				return _lam;
			}

			DynamicVector<Scalar> const& get_t() const
			{
				return _t;
			}

			//------------------------
			//
			// HPMPC memory interface
			//
			//------------------------
			double * x_data() { return _x.data(); }
			double * u_data() { return _u.data(); }
			double * pi_data() { return _pi.data(); }
			double * lam_data() { return _lam.data(); }
			double * t_data() { return _t.data(); }

		private:
			DynamicVector<Scalar> _x;
			DynamicVector<Scalar> _u;
			DynamicVector<Scalar> _pi;
			DynamicVector<Scalar> _lam;
			DynamicVector<Scalar> _t;
		};

		template <typename InputIterator>
		HPMPCSolution(InputIterator sz_first, InputIterator sz_last)
		{
			auto const N = std::distance(sz_first, sz_last);

			_stage.reserve(N);
			_x.reserve(N);
			_u.reserve(N);
			_pi.reserve(N);
			_lam.reserve(N);
			_t.reserve(N);

			for (; sz_first != sz_last; ++sz_first)
			{
				_stage.emplace_back(*sz_first, sz_first + 1 != sz_last ? sz_first[1].nx() : 0);
				Stage& st = _stage.back();

				_x.push_back(st.x_data());
				_u.push_back(st.u_data());
				_pi.push_back(st.pi_data());
				_lam.push_back(st.lam_data());
				_t.push_back(st.t_data());
			}
		}

		HPMPCSolution(HPMPCSolution const&) = delete;
		HPMPCSolution(HPMPCSolution&& rhs) = default;

		// ************************************************
		//                 HPMPC interface
		// ************************************************
		double * const * x_data() { return _x.data(); }
		double * const * u_data() { return _u.data(); }
		double * const * pi_data() { return _pi.data(); }
		double * const * lam_data() { return _lam.data(); }
		double * const * t_data() { return _t.data(); }
		double * inf_norm_res_data() { return _inf_norm_res.data(); }

		/// \brief Get number of iterations performed by the QP solver.
		unsigned getNumIter() const { return numIter_; }

		/// \brief Set number of iterations performed by the QP solver (called by the solver).
		void setNumIter(unsigned n) { numIter_ = n; }

	private:
		static Scalar constexpr sNaN = std::numeric_limits<Scalar>::signaling_NaN();

		std::vector<Stage> _stage;

		std::vector<double *> _x;
		std::vector<double *> _u;
		std::vector<double *> _pi;
		std::vector<double *> _lam;
		std::vector<double *> _t;
		StaticVector<double, 4> _inf_norm_res {sNaN};

		/// \brief Number of iterations performed by the QP solver.
		unsigned numIter_ = 0;
	};
}
