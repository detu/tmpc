/*
 * qpDUNESSolution.hpp
 *
 *  Created on: May 4, 2016
 *      Author: kotlyar
 */

#pragma once

//#include "../core/Trajectory.hpp"

#include <Eigen/Dense>

#include <vector>

namespace tmpc
{
	//
	// Provides solution interface for the HPMPC solver.
	//
	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	class HPMPCSolution //: public TrajectoryBase<HPMPCSolution<NX_, NU_>, NX_, NU_>
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

		StateVector const& get_x(std::size_t i) const
		{
			if (i >= nT() + 1)
				throw std::out_of_range("HPMPCSolution<>::get_x(): index is out of range");

			return i < nT() ? _stage[i]._x : _xEnd;
		}

		template <typename Matrix>
		void set_x(std::size_t i, Eigen::MatrixBase<Matrix> const& val)
		{
			stage(i)._x = val;
		}

		InputVector const& get_u(std::size_t i) const { return stage(i)._u; }
		StateVector const& get_xend() const	{ return _xEnd; }
		friend StateVector const& get_xend(HPMPCSolution const& s) { return s._xEnd; }

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

	private:
		struct StageData
		{
			StateVector _x;
			InputVector _u;
			StateVector _pi;
			LagrangeVector _lam;
			LagrangeVector _t;
		};

		std::vector<StageData> _stage;
		StateVector _xEnd;
		EndLagrangeVector _lamEnd;
		EndLagrangeVector _tEnd;

		std::vector<double *> _x;
		std::vector<double *> _u;
		std::vector<double *> _pi;
		std::vector<double *> _lam;
		std::vector<double *> _t;
		Eigen::Matrix<double, 4, 1> _inf_norm_res;

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
	};

	template<unsigned NX, unsigned NU, unsigned NC, unsigned NCT>
	typename HPMPCSolution<NX, NU, NC, NCT>::StateInputVector get_z(HPMPCSolution<NX, NU, NC, NCT> const& s, std::size_t i)
	{
		typename HPMPCSolution<NX, NU, NC, NCT>::StateInputVector z;
		z << s.get_x(i), s.get_u(i);
		return z;
	}
}
