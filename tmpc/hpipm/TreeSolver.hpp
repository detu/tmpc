#pragma once

#include <tmpc/hpipm/TreeOcpQpDim.hpp>
#include <tmpc/hpipm/TreeOcpQp.hpp>
#include <tmpc/hpipm/TreeOcpQpIpmArg.hpp>
#include <tmpc/hpipm/TreeOcpQpIpmWorkspace.hpp>
#include <tmpc/hpipm/TreeOcpQpSol.hpp>
#include <tmpc/hpipm/TreeOcpQpKkt.hpp>
#include <tmpc/hpipm/ConvertTreeOcpQp.hpp>

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/DynamicOcpSize.hpp>
#include <tmpc/ocp/OcpTree.hpp>
#include <tmpc/ocp/OcpSolution.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/property_map/PropertyMap.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/math/UnpaddedMatrix.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/Exception.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>
#include <memory>


namespace tmpc :: hpipm
{
	/**
	 * @brief Multistage QP solver using HPIPM
	 *
	 * @tparam <Real_> scalar floating-point type
	 */
	template <typename Real_>
	class TreeSolver
	{
	public:
		using Real = Real_;
		
		
		/**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <OcpSize Size>
		TreeSolver(Size const& size)
		:	graph_ {size.graph()}
		,	size_ {size}
		,	dim_ {size}
		,	solverArg_ {dim_, ROBUST}
		,	qp_ {dim_}
		,	sol_ {dim_}
		,	solverWorkspace_ {dim_, solverArg_}
		{
			auto const nv = num_vertices(graph_);
			auto const ne = num_edges(graph_);

			// Preallocate arrays holding data pointers.
			lbls_.reserve(nv);
			lbus_.reserve(nv);
			hidxbx_.reserve(nv);
			hidxbu_.reserve(nv);
			idxs_.reserve(nv);
	
			Zl_.reserve(nv);
			Zu_.reserve(nv);
			zl_.reserve(nv);
			zu_.reserve(nv);
			
			// Init vertex and edge data pointers.
			for (auto v : graph_.vertices())
				initVertex(v);
		}


		/**
		 * \brief Copying is not allowed.
		 */
		TreeSolver(TreeSolver const&) = delete;


		/// @brief Assignment prohibited.
		TreeSolver& operator= (TreeSolver const&) = delete;


		template <tmpc::OcpQp Qp, OcpSolution Solution>
		void operator()(Qp const& qp, Solution& sol)
		{
			convertQp(qp, qp_);
			solverWorkspace_.solve(qp_, sol_, solverArg_);
			convertSolution(sol);
		}


		template <tmpc::OcpQp Qp, OcpSolution Solution>
		void solveUnconstrained(Qp const& qp, Solution& sol)
		{
			convertQp(qp, qp_);
			solveUnconstrainedInternal();
			convertSolution(sol);
		}


		/// @brief Call HPIPM to solve an (unconstrained) QP which has already been copied to HPIPM structures
		/// by a previous call to convertProblem().
		void solveUnconstrainedInternal()
		{
			solverWorkspace_.fact_solve_kkt_unconstr(qp_, sol_, solverArg_);
		}


		/// @brief Copy QP solution from HPIPM data structures.
		template <OcpSolution Solution>
		void convertSolution(Solution& sol)
		{
			for (auto v : vertices(graph_))
			{
				sol_.get_u(v, sol.u(v).data());
				sol_.get_x(v, sol.x(v).data());

				// Slacks are not supported at the moment!
				// sol_->get_sl(v, get(sol.ls(), v).data());
				// sol_->get_su(v, get(sol.us(), v).data());

				blaze::DynamicVector<Real> lam_b(size_.nu(v) + size_.nx(v));

				// TODO: this works correctly only if all components of x and u are hard-constrained!
				sol_.get_lam_lb(v, lam_b.data());
				sol.lam_lu(v) = subvector(lam_b, 0, size_.nu(v));
				sol.lam_lx(v) = subvector(lam_b, size_.nu(v), size_.nx(v));
				
				sol_.get_lam_ub(v, lam_b.data());
				sol.lam_uu(v) = subvector(lam_b, 0, size_.nu(v));
				sol.lam_ux(v) = subvector(lam_b, size_.nu(v), size_.nx(v));

				sol_.get_lam_lg(v, sol.lam_ld(v).data());
				sol_.get_lam_ug(v, sol.lam_ud(v).data());
				
				// The following two functions should theretically be there, but they are not:
				// sol_->get_lam_sl(v, lam_ls_[v].data());
				// sol_->get_lam_su(v, lam_us_[v].data());
			}

			for (auto e : edges(graph_))
			{
				auto const v = source(e, graph_);
				sol_.get_pi(v, sol.pi(e).data());
			}
		}


		size_t maxIter() const noexcept 
		{ 
			return solverArg_->iter_max; 
		}

		
		/// \brief Get number of iterations performed by the QP solver.
		size_t numIter() const
		{
			return solverWorkspace_->get_iter();
		}


		void alphaMin(Real val)
		{
			if (val <= 0.)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Value must be positive"));

			solverArg_->set_alpha_min(val);
		}

		
		auto const& size() const noexcept
		{
			return size_;
		}


		auto const& graph() const
		{
			return graph_;
		}


	private:
		static auto constexpr SO = blaze::columnMajor;

		OcpTree graph_;
		DynamicOcpSize size_;

		// --------------------------------
		//
		// QP problem data
		//
		// --------------------------------
		std::vector<UnpaddedMatrix<Real, SO>> Zl_;
		std::vector<UnpaddedMatrix<Real, SO>> Zu_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> zl_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> zu_;

		// Lower bound of lower slack.
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lbls_;

		// Upper bound of upper slack.
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lbus_;

		// Hard constraints index for state
		std::vector<std::unique_ptr<int []>> hidxbx_;

		// Hard constraints index for input
		std::vector<std::unique_ptr<int []>> hidxbu_;

		// Soft constraints index
		std::vector<std::unique_ptr<int []>> idxs_;
		
		
		// --------------------------------
		//
		// QP solution data
		//
		// --------------------------------

		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> ls_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> us_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_ls_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_us_;

		// --------------------------------
		//
		// HPIPM solver data
		//
		// --------------------------------
		TreeOcpQpDim<Real> dim_;
		TreeOcpQpIpmArg<Real> solverArg_;
		TreeOcpQp<Real> qp_;
		TreeOcpQpSol<Real> sol_;
		TreeOcpQpIpmWorkspace<Real> solverWorkspace_;


		void initVertex(OcpVertex u)
		{
			{
				// Init idxbx to 0,1,2,... s.t. maximum necessary workspace size is returned.
				std::unique_ptr<int []> idxbx {new int[size_.nx(u)]};
				for (int i = 0; i < size_.nx(u); ++i)
					idxbx[i] = i;
				
				hidxbx_.emplace_back(std::move(idxbx));

				// Init idxbu to 0,1,2,... s.t. maximum necessary workspace size is returned.
				std::unique_ptr<int []> idxbu {new int [size_.nu(u)]};
				for (int i = 0; i < size_.nu(u); ++i)
					idxbu[i] = i;

				hidxbu_.emplace_back(std::move(idxbu));
			}
				
			{
				// Initialize idxs to 0,1,2, ... .
				std::unique_ptr<int []> idxs {new int[size_.ns(u)]};
				for (int i = 0; i < size_.ns(u); ++i)
					idxs[i] = i;

				idxs_.emplace_back(std::move(idxs));
			}

			Zl_.emplace_back(size_.ns(u), size_.ns(u));
			Zu_.emplace_back(size_.ns(u), size_.ns(u));
			zl_.emplace_back(size_.ns(u));
			zu_.emplace_back(size_.ns(u));				

			lbls_.emplace_back(size_.ns(u));
			lbus_.emplace_back(size_.ns(u));

			ls_.emplace_back(size_.ns(u));
			us_.emplace_back(size_.ns(u));
			
			lam_ls_.emplace_back(size_.ns(u));
			lam_us_.emplace_back(size_.ns(u));
		}
	};
}
