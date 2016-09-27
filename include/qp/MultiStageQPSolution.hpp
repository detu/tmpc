/*
 * MultiStageQPSolution.hpp
 *
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
	// Provides a generic solution interface for multistage QP solvers.
	//
	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	class MultiStageQPSolution
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
		typedef Eigen::Matrix<double, NC, 1> StageConstraintVector;
		typedef Eigen::Matrix<double, NC, 1> EndStageConstraintVector;
		typedef Eigen::Matrix<double, 2 * NC + 2 * (NX + NU), 1> LagrangeVector;
		typedef Eigen::Matrix<double, 2 * NCT + 2 * NX, 1> EndLagrangeVector;

		MultiStageQPSolution(size_type nt)
		:	_stage(nt)
		{
		}

		MultiStageQPSolution(MultiStageQPSolution const&) = delete;

		StateVector const& get_x(std::size_t i) const
		{
			if (i > nT())
				throw std::out_of_range("MultiStageQPSolution<>::get_x(): index is out of range");

			return i < nT() ? _stage[i]._x : _xEnd;
		}

		template <typename Matrix>
		void set_x(std::size_t i, Eigen::MatrixBase<Matrix> const& val)
		{
			if (i > nT())
				throw std::out_of_range("MultiStageQPSolution<>::set_x(): index is out of range");

			(i < nT() ? _stage[i]._x : _xEnd) = val;
		}

		template <typename Matrix>
		void set_u(std::size_t i, Eigen::MatrixBase<Matrix> const& val)
		{
			stage(i)._u = val;
		}

		InputVector const& get_u(std::size_t i) const { return stage(i)._u; }

		size_type constexpr nX() { return NX; }
		size_type constexpr nU() { return NU; }
		size_type nT() const noexcept { return _stage.size(); }

		StateVector const& get_pi(std::size_t i) const { return stage(i)._pi; }

		decltype(auto) get_lam_u_min(std::size_t i) const { return signaling_nan<InputVector>(); }
		decltype(auto) get_lam_u_max(std::size_t i) const { return signaling_nan<InputVector>(); }
		decltype(auto) get_lam_x_min(std::size_t i) const {	return signaling_nan<StateVector>(); }
		decltype(auto) get_lam_x_max(std::size_t i) const {	return signaling_nan<StateVector>(); }

		decltype(auto) get_lam_d_min(std::size_t i) const { return signaling_nan<StageConstraintVector>(); }
		decltype(auto) get_lam_d_max(std::size_t i) const { return signaling_nan<StageConstraintVector>(); }

		decltype(auto) get_lam_d_end_min() const { return signaling_nan<EndStageConstraintVector>(); }
		decltype(auto) get_lam_d_end_max() const { return signaling_nan<EndStageConstraintVector>(); }

	private:
		struct StageData
		{
			StateVector    _x   = signaling_nan<StateVector   >();
			InputVector    _u   = signaling_nan<InputVector   >();
			StateVector    _pi  = signaling_nan<StateVector   >();
			LagrangeVector _lam = signaling_nan<LagrangeVector>();
			LagrangeVector _t   = signaling_nan<LagrangeVector>();
		};

		std::vector<StageData> _stage;
		StateVector       _xEnd   = signaling_nan<StateVector      >();
		EndLagrangeVector _lamEnd = signaling_nan<EndLagrangeVector>();
		EndLagrangeVector _tEnd   = signaling_nan<EndLagrangeVector>();

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
}
