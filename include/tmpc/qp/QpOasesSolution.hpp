/*
 * MultiStageQPSolution.hpp
 *
 *      Author: kotlyar
 */

#pragma once

#include "QpSize.hpp"

#include <tmpc/Matrix.hpp>

#include <vector>

namespace tmpc
{
	/**
	 * \brief QP solution with memory layout matching qpOASES
	 */
	class QpOasesSolution
	{
	public:
		typedef std::size_t size_type;
		typedef CustomVector<double, unaligned, unpadded> VectorMap;

		template <typename InputIt>
		QpOasesSolution(InputIt sz_begin, InputIt sz_end)
		:	primalSolution_(numVariables(sz_begin, sz_end))
		,	dualSolution_(numVariables(sz_begin, sz_end) + numEqualities(sz_begin, sz_end) + numInequalities(sz_begin, sz_end))
		,	size_(sz_begin, sz_end)
		{
			assert(stage_.empty());
			stage_.reserve(size_.size());

			double * x = primalSolution_.data();
			double * lambda_x = dualSolution_.data();
			double * lambda = dualSolution_.data() + primalSolution_.size();

			for (auto s = sz_begin; s != sz_end; ++s)
			{
				auto const s_next = s + 1;
				auto const nx_next = s_next != sz_end ? s_next->nx() : 0;
				stage_.emplace_back(*s, nx_next, x, lambda_x, lambda);

				x += s->nx() + s->nu();
				lambda_x += s->nx() + s->nu();
				lambda += s->nc() + nx_next;
			}
		}

		QpOasesSolution(std::vector<QpSize> const& sz);
		QpOasesSolution(QpOasesSolution const&) = delete;
		QpOasesSolution(QpOasesSolution &&) = default;

		//---------------------------
		// Access interface
		//---------------------------

		class Stage
		{
		public:
			Stage(QpSize const& sz, std::size_t nx_plus, double * x, double * lam_x, double * lam);

			VectorMap const& get_x() const
			{
				return x_;
			}

			VectorMap const& get_u() const
			{
				return u_;
			}

			VectorMap const& get_pi() const
			{
				return pi_;
			}

			VectorMap const& get_lam_u() const
			{
				return lamU_;
			}

			VectorMap const& get_lam_x() const
			{
				return lamX_;
			}

			VectorMap const& get_lam_d() const
			{
				return lam_;
			}

		private:
			VectorMap x_;
			VectorMap u_;
			VectorMap lamX_;
			VectorMap lamU_;
			VectorMap lam_;
			VectorMap pi_;
		};

		Stage const& operator[](std::size_t i) const
		{
			return stage_.at(i);
		}

		Stage& operator[](std::size_t i)
		{
			return stage_.at(i);
		}

		typedef std::vector<Stage>::iterator iterator;
		typedef std::vector<Stage>::const_iterator const_iterator;
		typedef std::vector<Stage>::reference reference;
		typedef std::vector<Stage>::const_reference const_reference;

		iterator begin()
		{
			return stage_.begin();
		}

		iterator end()
		{
			return stage_.end();
		}

		const_iterator begin() const
		{
			return stage_.begin();
		}

		const_iterator end() const
		{
			return stage_.end();
		}

		reference front()
		{
			return stage_.front();
		}

		reference back()
		{
			return stage_.back();
		}

		const_reference front() const
		{
			return stage_.front();
		}

		const_reference back() const
		{
			return stage_.back();
		}

		//---------------------------
		// Obsolete access interface
		//---------------------------

		VectorMap const& get_x(std::size_t i) const
		{
			return stage(i).get_x();
		}

		VectorMap const& get_u(std::size_t i) const
		{
			return stage(i).get_u();
		}

		VectorMap const& get_pi(std::size_t i) const
		{
			return stage(i).get_pi();

		}

		decltype(auto) get_lam_u_min(std::size_t i) const
		{
			return forEach(stage(i).get_lam_u(), [] (double d) { return max(d, 0.); });
		}

		decltype(auto) get_lam_u_max(std::size_t i) const
		{
			return -forEach(stage(i).get_lam_u(), [] (double d) { return min(d, 0.); });
		}

		decltype(auto) get_lam_x_min(std::size_t i) const
		{
			return forEach(stage(i).get_lam_x(), [] (double d) { return max(d, 0.); });
		}

		decltype(auto) get_lam_x_max(std::size_t i) const
		{
			return -forEach(stage(i).get_lam_x(), [] (double d) { return max(d, 0.); });
		}

		decltype(auto) get_lam_d_min(std::size_t i) const
		{
			return forEach(stage(i).get_lam_d(), [] (double d) { return max(d, 0.); });
		}

		decltype(auto) get_lam_d_max(std::size_t i) const
		{
			return -forEach(stage(i).get_lam_d(), [] (double d) { return max(d, 0.); });
		}

		decltype(auto) get_lam_d_end_min() const
		{
			return forEach(stage_.back().get_lam_d(), [] (double d) { return max(d, 0.); });
		}

		decltype(auto) get_lam_d_end_max() const
		{
			return -forEach(stage_.back().get_lam_d(), [] (double d) { return max(d, 0.); });
		}

		/// \brief Get number of iterations performed by the QP solver.
		unsigned getNumIter() const { return numIter_; }

		/// \brief Set number of iterations performed by the QP solver (called by the solver).
		void setNumIter(unsigned n) { numIter_ = n; }

		// ****************** qpOASES memory interface *******************

		double * primalSolutionData()
		{
			return primalSolution_.data();
		}

		double * dualSolutionData()
		{
			return dualSolution_.data();
		}

		std::vector<QpSize> const& size() const
		{
			return size_;
		}

	private:
		Stage const& stage(std::size_t i) const
		{
			return stage_.at(i);
		}

		Stage& stage(std::size_t i)
		{
			return stage_.at(i);
		}

		DynamicVector<double> primalSolution_;
		DynamicVector<double> dualSolution_;

		std::vector<QpSize> size_;
		std::vector<Stage> stage_;

		/// \brief Number of iterations performed by the QP solver.
		unsigned numIter_ = 0;
	};
}
