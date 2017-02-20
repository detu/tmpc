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
	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	class HPMPCSolution
	{
	public:
		typedef std::size_t size_type;

		static size_type const NX = NX_;
		static size_type const NU = NU_;
		static size_type const NZ = NX + NU;
		static size_type const NC = NC_;
		static size_type const NCT = NCT_;

		typedef StaticVector<double, NX> StateVector;
		typedef StaticVector<double, NU> InputVector;
		typedef StaticVector<double, NZ> StateInputVector;
		typedef StaticVector<double, 2 * NC + 2 * (NX + NU)> LagrangeVector;
		typedef StaticVector<double, 2 * NCT + 2 * NX> EndLagrangeVector;

		HPMPCSolution(size_type nt)
		:	_stage(nt)
		,	_x    (nt + 1)
		,	_u    (nt    )
		,	_pi   (nt    )
		,	_lam  (nt + 1)
		,	_t    (nt + 1)
		{
			for (std::size_t i = 0; i < nt; ++i)
			{
				_x[i]   = _stage[i]._x  .data();
				_u[i]   = _stage[i]._u  .data();
				_pi[i]  = _stage[i]._pi .data();
				_lam[i] = _stage[i]._lam.data();
				_t[i]   = _stage[i]._t  .data();
			}

			_x  .back() = _xEnd  .data();
			_lam.back() = _lamEnd.data();
			_t  .back() = _tEnd  .data();
		}

		HPMPCSolution(HPMPCSolution const&) = delete;

		HPMPCSolution(HPMPCSolution&& rhs)
		:	_stage(std::move(rhs._stage))
		,	_xEnd(std::move(rhs._xEnd))
		,	_lamEnd(std::move(rhs._lamEnd))
		,	_tEnd(std::move(rhs._tEnd))
		,	_x    (std::move(rhs._x))
		,	_u    (std::move(rhs._u))
		,	_pi   (std::move(rhs._pi))
		,	_lam  (std::move(rhs._lam))
		,	_t    (std::move(rhs._t))
		,	_inf_norm_res(std::move(rhs._inf_norm_res))
		,	numIter_(std::move(rhs.numIter_))
		{
			_x  .back() = _xEnd  .data();
			_lam.back() = _lamEnd.data();
			_t  .back() = _tEnd  .data();
		}

		StateVector const& get_x(std::size_t i) const
		{
			if (i > nT())
				throw std::out_of_range("HPMPCSolution<>::get_x(): index is out of range");

			return i < nT() ? _stage[i]._x : _xEnd;
		}

		template <typename Matrix>
		void set_x(std::size_t i, Matrix const& val)
		{
			if (i > nT())
				throw std::out_of_range("HPMPCSolution<>::set_x(): index is out of range");

			(i < nT() ? _stage[i]._x : _xEnd) = val;
		}

		InputVector const& get_u(std::size_t i) const { return stage(i)._u; }

		StateVector const& get_pi(std::size_t i) const { return stage(i)._pi; }

		decltype(auto) get_lam_u_min(std::size_t i) const
		{
			return subvector(stage(i)._lam, 0., NU);
		}

		decltype(auto) get_lam_u_max(std::size_t i) const
		{
			return subvector(stage(i)._lam, NU + NX + NC, NU);
		}

		CustomVector<double, unaligned, unpadded> get_lam_x_min(std::size_t i) const
		{
			if (i > nT())
				throw std::out_of_range("HPMPCSolution::get_lam_x_min(): index is out of range");

			return CustomVector<double, unaligned, unpadded>(_lam[i] + NU, NX);
		}

		CustomVector<double, unaligned, unpadded> get_lam_x_max(std::size_t i) const
		{
			if (i > nT())
				throw std::out_of_range("HPMPCSolution::get_lam_x_max(): index is out of range");

			return CustomVector<double, unaligned, unpadded>(_lam[i] + NU + NX + NC + NU, NX);
		}

		Subvector<LagrangeVector> get_lam_d_min(std::size_t i) const
		{
			return subvector(stage(i)._lam, NU + NX, NC);
		}

		Subvector<LagrangeVector> get_lam_d_max(std::size_t i) const
		{
			return subvector(stage(i)._lam, NU + NX + NC + NU + NX, NC);
		}

		Subvector<EndLagrangeVector> get_lam_d_end_min() const
		{
			return subvector(_lamEnd, NX, NC);
		}

		Subvector<EndLagrangeVector> get_lam_d_end_max() const
		{
			return subvector(_lamEnd, NX + NCT + NX, NC);
		}

		size_type const nX() const noexcept { return NX; }
		size_type const nU() const noexcept { return NU; }
		size_type const nT() const { return _stage.size(); }

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
		struct StageData
		{
			// Initialize all numeric data to NaN so that if an uninitialized object
			// by mistake used in calculations is easier to detect.
			StateVector      _x = sNaN;
			InputVector      _u = sNaN;
			StateVector     _pi = sNaN;
			LagrangeVector _lam = sNaN;
			LagrangeVector   _t = sNaN;
		};

		StageData& stage(size_type i)
		{
			assert(i < nT());
			return _stage[i];
		}

		StageData const& stage(size_type i) const
		{
			assert(i < nT());
			return _stage[i];
		}

		static double constexpr sNaN = std::numeric_limits<double>::signaling_NaN();

		std::vector<StageData> _stage;

		// Initialize all numeric data to NaN so that if an uninitialized object
		// by mistake used in calculations is easier to detect.
		StateVector       _xEnd   = sNaN;
		EndLagrangeVector _lamEnd = sNaN;
		EndLagrangeVector _tEnd   = sNaN;

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
