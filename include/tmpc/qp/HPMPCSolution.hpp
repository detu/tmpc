/*
 * qpDUNESSolution.hpp
 *
 *  Created on: May 4, 2016
 *      Author: kotlyar
 */

#pragma once

#include "qp.hpp"
#include "../core/matrix.hpp"

#include <Eigen/Dense>

#include <vector>

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

		typedef Eigen::Matrix<double, NX, 1> StateVector;
		typedef Eigen::Matrix<double, NU, 1> InputVector;
		typedef Eigen::Matrix<double, NZ, 1> StateInputVector;
		typedef Eigen::Matrix<double, 2 * NC + 2 * (NX + NU), 1> LagrangeVector;
		typedef Eigen::Matrix<double, 2 * NCT + 2 * NX, 1> EndLagrangeVector;

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
		void set_x(std::size_t i, Eigen::MatrixBase<Matrix> const& val)
		{
			if (i > nT())
				throw std::out_of_range("HPMPCSolution<>::set_x(): index is out of range");

			(i < nT() ? _stage[i]._x : _xEnd) = val;
		}

		InputVector const& get_u(std::size_t i) const { return stage(i)._u; }

		StateVector const& get_pi(std::size_t i) const { return stage(i)._pi; }

		decltype(auto) get_lam_u_min(std::size_t i) const { return middle_rows<NU>(stage(i)._lam, 0); }
		decltype(auto) get_lam_u_max(std::size_t i) const { return middle_rows<NU>(stage(i)._lam, NU + NX + NC); }

		Eigen::Map<StateVector> get_lam_x_min(std::size_t i) const
		{
			if (i > nT())
				throw std::out_of_range("HPMPCSolution::get_lam_x_min(): index is out of range");

			return Eigen::Map<StateVector>(_lam[i] + NU);
		}

		Eigen::Map<StateVector> get_lam_x_max(std::size_t i) const
		{
			if (i > nT())
				throw std::out_of_range("HPMPCSolution::get_lam_x_max(): index is out of range");

			return Eigen::Map<StateVector>(_lam[i] + NU + NX + NC + NU);
		}

		decltype(auto) get_lam_d_min(std::size_t i) const { return middle_rows<NC>(stage(i)._lam, NU + NX); }
		decltype(auto) get_lam_d_max(std::size_t i) const { return middle_rows<NC>(stage(i)._lam, NU + NX + NC + NU + NX); }

		decltype(auto) get_lam_d_end_min() const { return middle_rows<NC>(_lamEnd, NX); }
		decltype(auto) get_lam_d_end_max() const { return middle_rows<NC>(_lamEnd, NX + NCT + NX); }

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
			StateVector      _x = signaling_nan<StateVector   >();
			InputVector      _u = signaling_nan<InputVector   >();
			StateVector     _pi = signaling_nan<StateVector   >();
			LagrangeVector _lam = signaling_nan<LagrangeVector>();
			LagrangeVector   _t = signaling_nan<LagrangeVector>();
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

		std::vector<StageData> _stage;

		// Initialize all numeric data to NaN so that if an uninitialized object
		// by mistake used in calculations is easier to detect.
		StateVector       _xEnd   = signaling_nan<StateVector      >();
		EndLagrangeVector _lamEnd = signaling_nan<EndLagrangeVector>();
		EndLagrangeVector _tEnd   = signaling_nan<EndLagrangeVector>();

		std::vector<double *> _x;
		std::vector<double *> _u;
		std::vector<double *> _pi;
		std::vector<double *> _lam;
		std::vector<double *> _t;
		Eigen::Matrix<double, 4, 1> _inf_norm_res = signaling_nan<Eigen::Matrix<double, 4, 1>>();

		/// \brief Number of iterations performed by the QP solver.
		unsigned numIter_ = 0;
	};
}
