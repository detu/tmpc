/*
 * MultiStageQPSolution.hpp
 *
 *      Author: kotlyar
 */

#pragma once

#include "QpSize.hpp"

#include <Eigen/Dense>

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
		typedef Eigen::Map<Eigen::VectorXd> VectorMap;

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

		VectorMap const& get_x(std::size_t i) const
		{
			return stage(i).x_;
		}

		template <typename Matrix>
		void set_x(std::size_t i, Eigen::MatrixBase<Matrix> const& val)
		{
			stage(i).x_ = val;
		}

		VectorMap const& get_u(std::size_t i) const
		{
			return stage(i).u_;
		}

		template <typename Matrix>
		void set_u(std::size_t i, Eigen::MatrixBase<Matrix> const& val)
		{
			stage(i).u_ = val;
		}

		VectorMap const& get_pi(std::size_t i) const
		{
			return stage(i).pi_;
		}

		decltype(auto) get_lam_u_min(std::size_t i) const
		{
			return stage(i).lamU_.array().max(0.).matrix();
		}

		decltype(auto) get_lam_u_max(std::size_t i) const
		{
			return -stage(i).lamU_.array().min(0.).matrix();
		}

		decltype(auto) get_lam_x_min(std::size_t i) const
		{
			return stage(i).lamX_.array().max(0.).matrix();
		}

		decltype(auto) get_lam_x_max(std::size_t i) const
		{
			return -stage(i).lamX_.array().min(0.).matrix();
		}

		decltype(auto) get_lam_d_min(std::size_t i) const
		{
			return stage(i).lam_.array().max(0.).matrix();
		}

		decltype(auto) get_lam_d_max(std::size_t i) const
		{
			return -stage(i).lam_.array().min(0.).matrix();
		}

		decltype(auto) get_lam_d_end_min() const
		{
			return stage_.back().lam_.array().max(0.).matrix();
		}

		decltype(auto) get_lam_d_end_max() const
		{
			return -stage_.back().lam_.array().min(0.).matrix();
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
		struct Stage
		{
			Stage(QpSize const& sz, std::size_t nx_plus, double * x, double * lam_x, double * lam);

			VectorMap x_;
			VectorMap u_;
			VectorMap lamX_;
			VectorMap lamU_;
			VectorMap lam_;
			VectorMap pi_;
		};

		Stage const& stage(std::size_t i) const
		{
			return stage_.at(i);
		}

		Stage& stage(std::size_t i)
		{
			return stage_.at(i);
		}

		Eigen::VectorXd primalSolution_;
		Eigen::VectorXd dualSolution_;

		std::vector<QpSize> size_;
		std::vector<Stage> stage_;

		/// \brief Number of iterations performed by the QP solver.
		unsigned numIter_ = 0;
	};
}
