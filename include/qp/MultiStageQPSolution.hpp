/*
 * MultiStageQPSolution.hpp
 *
 *      Author: kotlyar
 */

#pragma once

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
		typedef Eigen::Matrix<double, 2 * NC + 2 * (NX + NU), 1> LagrangeVector;
		typedef Eigen::Matrix<double, 2 * NCT + 2 * NX, 1> EndLagrangeVector;

		MultiStageQPSolution(size_type nt)
		:	_stage(nt)
		{
		}

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

		size_type const nX() const noexcept { return NX; }
		size_type const nU() const noexcept { return NU; }
		size_type const nT() const { return _stage.size(); }

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
