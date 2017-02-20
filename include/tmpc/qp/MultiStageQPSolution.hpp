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
	template <unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	class MultiStageQPSolution
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
		typedef StaticVector<double, NC> StageConstraintVector;
		typedef StaticVector<double, NC> EndStageConstraintVector;
		typedef StaticVector<double, 2 * NC + 2 * (NX + NU)> LagrangeVector;
		typedef StaticVector<double, 2 * NCT + 2 * NX> EndLagrangeVector;

		MultiStageQPSolution(size_type nt)
		:	_stage(nt)
		{
		}

		MultiStageQPSolution(MultiStageQPSolution const&) = delete;
		MultiStageQPSolution(MultiStageQPSolution&&) = default;

		StateVector const& get_x(std::size_t i) const
		{
			if (i > nT())
				throw std::out_of_range("MultiStageQPSolution<>::get_x(): index is out of range");

			return i < nT() ? _stage[i]._x : _xEnd;
		}

		template <typename Matrix>
		void set_x(std::size_t i, Matrix const& val)
		{
			if (i > nT())
				throw std::out_of_range("MultiStageQPSolution<>::set_x(): index is out of range");

			(i < nT() ? _stage[i]._x : _xEnd) = val;
		}

		template <typename Matrix>
		void set_u(std::size_t i, Matrix const& val)
		{
			stage(i)._u = val;
		}

		InputVector const& get_u(std::size_t i) const { return stage(i)._u; }

		size_type constexpr nX() { return NX; }
		size_type constexpr nU() { return NU; }
		size_type nT() const noexcept { return _stage.size(); }

		StateVector const& get_pi(std::size_t i) const { return stage(i)._pi; }

		decltype(auto) get_lam_u_min(std::size_t i) const { return InputVector(sNaN); }
		decltype(auto) get_lam_u_max(std::size_t i) const { return InputVector(sNaN); }
		decltype(auto) get_lam_x_min(std::size_t i) const {	return StateVector(sNaN); }
		decltype(auto) get_lam_x_max(std::size_t i) const {	return StateVector(sNaN); }

		decltype(auto) get_lam_d_min(std::size_t i) const { return StageConstraintVector(sNaN); }
		decltype(auto) get_lam_d_max(std::size_t i) const { return StageConstraintVector(sNaN); }

		decltype(auto) get_lam_d_end_min() const { return EndStageConstraintVector(sNaN); }
		decltype(auto) get_lam_d_end_max() const { return EndStageConstraintVector(sNaN); }

		/// \brief Get number of iterations performed by the QP solver.
		unsigned getNumIter() const { return numIter_; }

		/// \brief Set number of iterations performed by the QP solver (called by the solver).
		void setNumIter(unsigned n) { numIter_ = n; }

	private:
		static constexpr double sNaN = std::numeric_limits<double>::signaling_NaN();

		struct StageData
		{
			StateVector    _x   { sNaN };
			InputVector    _u   { sNaN };
			StateVector    _pi  { sNaN };
			LagrangeVector _lam { sNaN };
			LagrangeVector _t   { sNaN };
		};

		std::vector<StageData> _stage;
		StateVector       _xEnd   { sNaN };
		EndLagrangeVector _lamEnd { sNaN };
		EndLagrangeVector _tEnd   { sNaN };

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

		/// \brief Number of iterations performed by the QP solver.
		unsigned numIter_ = 0;
	};
}
