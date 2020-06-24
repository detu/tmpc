#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/math/UnpaddedMatrix.hpp>
#include <tmpc/qp/hpipm/OcpQpDim.hpp>
#include <tmpc/qp/hpipm/OcpQp.hpp>
#include <tmpc/qp/hpipm/OcpQpIpmArg.hpp>
#include <tmpc/qp/hpipm/OcpQpIpmWs.hpp>
#include <tmpc/qp/hpipm/OcpQpSol.hpp>
#include <tmpc/qp/hpipm/OcpQpKkt.hpp>
#include <tmpc/Traits.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/Exception.hpp>

#include <tmpc/property_map/VectorPtrPropertyMap.hpp>
#include <tmpc/property_map/MatrixPtrPropertyMap.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>
#include <memory>


namespace tmpc
{
	/**
	 * \brief Multistage QP solver using HPIPM
	 *
	 * \tparam <Kernel_> matrix kernel type
	 */
	template <typename Real>
	class HpipmWorkspace
	{
	public:
		static auto constexpr SO = blaze::columnMajor;
		

		std::string impl_solverName() const
		{
			return "HPIPM";
		}

	
		/**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <typename SizeMap>
		HpipmWorkspace(OcpGraph const& graph, SizeMap size, unsigned max_iter = 100)
		:	graph_ {graph}
		,	size_(num_vertices(graph))
		,	ocpQpDim_(makeOcpQpDim(graph, size))
		,	solverArg_ {new hpipm::OcpQpIpmArg<Real> {*ocpQpDim_, ROBUST}}
		,	qp_ {new hpipm::OcpQp<Real> {*ocpQpDim_}}
		,	sol_ {new hpipm::OcpQpSol<Real> {*ocpQpDim_}}
		,	solverWorkspace_ {new hpipm::OcpQpIpmWs<Real> {*ocpQpDim_, *solverArg_}}
		{
			copyProperty(size, iterator_property_map(size_.begin(), get(graph::vertex_index, graph_)), 
				graph::vertices(graph_));
			
			auto const nv = num_vertices(graph_);
			auto const ne = num_edges(graph_);

			if (nv < 2)
				throw std::invalid_argument("HPIPM needs an at least 2-stages problem");

			if (ne + 1 != num_vertices(graph_))
				throw std::invalid_argument("Number of edges in HPIPM size graph must be 1 less than the number of vertices");

			// solverArg_->alpha_min = 1e-8;
			// solverArg_->res_g_max = 1e-8;
			// solverArg_->res_b_max = 1e-8;
			// solverArg_->res_d_max = 1e-12;
			// solverArg_->res_m_max = 1e-12;
			// solverArg_->mu0 = 2.0;
			// solverArg_->iter_max = max_iter;
			// solverArg_->stat_max = 100;
			// solverArg_->pred_corr = 1;
			
			// Preallocate arrays holding QP vertex data pointers.
			lbls_.reserve(nv);
			lbus_.reserve(nv);
			hidxbx_.reserve(nv);
			hidxbu_.reserve(nv);
			idxs_.reserve(nv);
	
			Q_ .reserve(nv);
			S_ .reserve(nv);
			R_ .reserve(nv);
			q_ .reserve(nv);
			r_ .reserve(nv);
			Zl_.reserve(nv);
			Zu_.reserve(nv);
			zl_.reserve(nv);
			zu_.reserve(nv);
			lx_.reserve(nv);
			ux_.reserve(nv);
			lu_.reserve(nv);
			uu_.reserve(nv);
			C_ .reserve(nv);
			D_ .reserve(nv);
			lg_.reserve(nv);
			ug_.reserve(nv);
			
			x_.reserve(nv);
			u_.reserve(nv);
			lam_lb_.reserve(nv);
			lam_ub_.reserve(nv);
			lam_lg_.reserve(nv);
			lam_ug_.reserve(nv);

			// Preallocate arrays holding QP edge data pointers.
			A_ .reserve(ne);
			B_ .reserve(ne);
			b_ .reserve(ne);
			pi_.reserve(ne);

			// Traverse the graph and init vertex and edge data pointers.
			breadth_first_search(graph_, vertex(0, graph_), visitor(InitVisitor(*this)));
		}


		/**
		 * \brief Copy constructor
		 *
		 * Copying is not allowed.
		 */
		HpipmWorkspace(HpipmWorkspace const&) = delete;

		/**
		 * \brief Move constructor
		 *
		 * Move-construction is ok.
		 */
		HpipmWorkspace(HpipmWorkspace&& rhs) = default;

		HpipmWorkspace& operator= (HpipmWorkspace const&) = delete;
		HpipmWorkspace& operator= (HpipmWorkspace &&) = delete;

		void solve()
		{
			if (num_vertices(graph_) > 1)
			{
				convertProblem();
				hpipm::ocp_qp_ipm_solve(*qp_, *sol_, *solverArg_, *solverWorkspace_);
				convertSolution();
			}
		}


		void solveUnconstrained()
		{
			convertProblem();
			solveUnconstrainedInternal();
			convertSolution();
		}


		/// @brief Call HPIPM to solve an (unconstrained) QP which has already been copied to HPIPM structures
		/// by a previous call to convertProblem().
		void solveUnconstrainedInternal()
		{
			hpipm::fact_solve_kkt_unconstr_ocp_qp(*qp_, *sol_, *solverArg_, *solverWorkspace_);
		}


		/// @brief Copy QP to HPIPM data structures.
		void convertProblem()
		{
			for (auto v : graph::vertices(graph_))
			{
				qp_->set_Q(v, Q_[v].data());
				qp_->set_R(v, R_[v].data());
				qp_->set_S(v, S_[v].data());
				qp_->set_q(v, q_[v].data());
				qp_->set_r(v, r_[v].data());
				qp_->set_lbx(v, data(lx_[v]));
				qp_->set_ubx(v, data(ux_[v]));
				qp_->set_idxbx(v, hidxbx_[v].get());
				qp_->set_lbu(v, data(lu_[v]));
				qp_->set_ubu(v, data(uu_[v]));
				qp_->set_idxbu(v, hidxbu_[v].get());
				qp_->set_C(v, C_[v].data());
				qp_->set_D(v, D_[v].data());
				qp_->set_lg(v, lg_[v].data());
				qp_->set_ug(v, ug_[v].data());
				qp_->set_Zl(v, Zl_[v].data());
				qp_->set_Zu(v, Zu_[v].data());
				qp_->set_zl(v, zl_[v].data());
				qp_->set_zu(v, zu_[v].data());
				qp_->set_idxs(v, idxs_[v].get());
				qp_->set_lls(v, lbls_[v].data());
				qp_->set_lus(v, lbus_[v].data());
			}

			for (auto e : graph::edges(graph_))
			{
				auto const v = source(e, graph_);

				qp_->set_A(v, A_[v].data());
				qp_->set_B(v, B_[v].data());
				qp_->set_b(v, b_[v].data());
			}
		}


		/// @brief Copy QP solution from HPIPM data structures.
		void convertSolution()
		{
			for (auto v : graph::vertices(graph_))
			{
				sol_->get_u(v, u_[v].data());
				sol_->get_x(v, x_[v].data());
				sol_->get_sl(v, ls_[v].data());
				sol_->get_su(v, us_[v].data());
				sol_->get_lam_lb(v, lam_lb_[v].data());
				sol_->get_lam_ub(v, lam_ub_[v].data());
				sol_->get_lam_lg(v, lam_lg_[v].data());
				sol_->get_lam_ug(v, lam_ug_[v].data());
				sol_->get_lam_ub(v, lam_ub_[v].data());
				
				// The following two functions should theretically be there, but they are not:
				// sol_->get_lam_sl(v, lam_ls_[v].data());
				// sol_->get_lam_su(v, lam_us_[v].data());
			}

			for (auto e : graph::edges(graph_))
			{
				auto const v = source(e, graph_);

				sol_->get_pi(v, pi_[v].data());
			}
		}


		size_t maxIter() const noexcept 
		{ 
			return solverArg_->iter_max; 
		}

		
		void maxIter(size_t val) 
		{ 
			if (val != maxIter())
			{
				solverArg_->set_iter_max(val);
				recreateSolverWorkspace();
			}
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

		
		// void resGMax(Real val)
		// {
		// 	if (val <= 0.)
		// 		throw std::invalid_argument("HpipmWorkspace::resGMax(): value must be positive");

		// 	solverArg_->res_g_max = val;
		// }

		
		// void resBMax(Real val)
		// {
		// 	if (val <= 0.)
		// 		throw std::invalid_argument("HpipmWorkspace::resBMax(): value must be positive");
			
		// 	solverArg_->res_b_max = val;
		// }

		
		// void resDMax(Real val)
		// {
		// 	if (val <= 0.)
		// 		throw std::invalid_argument("HpipmWorkspace::resDMax(): value must be positive");
			
		// 	solverArg_->res_d_max = val;
		// }

		
		// void resMMax(Real val)
		// {
		// 	if (val <= 0.)
		// 		throw std::invalid_argument("HpipmWorkspace::resMMax(): value must be positive");
		
		// 	solverArg_->res_m_max = val;
		// }


		// void mu0(Real val)
		// {
		// 	if (val <= 0.)
		// 		throw std::invalid_argument("HpipmWorkspace::mu0(): value must be positive");
		
		// 	solverArg_->mu0 = val;
		// }
		
		
		auto size() const
		{
			return iterator_property_map(size_.begin(), get(graph::vertex_index, graph_));
		}


		auto const& graph() const
		{
			return graph_;
		}


		auto Q()
		{
			return make_iterator_property_map(Q_.begin(), get(graph::vertex_index, graph_));
		}


		auto Q() const
		{
			return make_iterator_property_map(Q_.begin(), get(graph::vertex_index, graph_));
		}


		auto R()
		{
			return make_iterator_property_map(R_.begin(), get(graph::vertex_index, graph_));
		}


		auto R() const
		{
			return make_iterator_property_map(R_.begin(), get(graph::vertex_index, graph_));
		}


		auto S()
		{
			return make_iterator_property_map(S_.begin(), get(graph::vertex_index, graph_));
		}


		auto S() const
		{
			return make_iterator_property_map(S_.begin(), get(graph::vertex_index, graph_));
		}


		auto q()
		{
			return make_iterator_property_map(q_.begin(), get(graph::vertex_index, graph_));
		}


		auto q() const
		{
			return make_iterator_property_map(q_.begin(), get(graph::vertex_index, graph_));
		}


		auto r()
		{
			return make_iterator_property_map(r_.begin(), get(graph::vertex_index, graph_));
		}


		auto r() const
		{
			return make_iterator_property_map(r_.begin(), get(graph::vertex_index, graph_));
		}


		auto lx()
		{
			// TODO: check size of put()
			return make_iterator_property_map(lx_.begin(), get(graph::vertex_index, graph_));
		}


		auto lx() const
		{
			// TODO: check size of put()
			return make_iterator_property_map(lx_.begin(), get(graph::vertex_index, graph_));
		}


		auto ux()
		{
			// TODO: check size of put()
			return make_iterator_property_map(ux_.begin(), get(graph::vertex_index, graph_));
		}


		auto ux() const
		{
			// TODO: check size of put()
			return make_iterator_property_map(ux_.begin(), get(graph::vertex_index, graph_));
		}


		auto lu()
		{
			// TODO: check size of put()
			return make_iterator_property_map(lu_.begin(), get(graph::vertex_index, graph_));
		}


		auto lu() const
		{
			// TODO: check size of put()
			return make_iterator_property_map(lu_.begin(), get(graph::vertex_index, graph_));
		}


		auto uu()
		{
			// TODO: check size of put()
			return make_iterator_property_map(uu_.begin(), get(graph::vertex_index, graph_));
		}


		auto uu() const
		{
			// TODO: check size of put()
			return make_iterator_property_map(uu_.begin(), get(graph::vertex_index, graph_));
		}


		auto C()
		{
			return make_iterator_property_map(C_.begin(), get(graph::vertex_index, graph_));
		}


		auto C() const
		{
			return make_iterator_property_map(C_.begin(), get(graph::vertex_index, graph_));
		}


		auto D()
		{
			return make_iterator_property_map(D_.begin(), get(graph::vertex_index, graph_));
		}


		auto D() const
		{
			return make_iterator_property_map(D_.begin(), get(graph::vertex_index, graph_));
		}


		auto ld()
		{
			return make_iterator_property_map(lg_.begin(), get(graph::vertex_index, graph_));
		}


		auto ld() const
		{
			return make_iterator_property_map(lg_.begin(), get(graph::vertex_index, graph_));
		}


		auto ud()
		{
			return make_iterator_property_map(ug_.begin(), get(graph::vertex_index, graph_));
		}


		auto ud() const
		{
			return make_iterator_property_map(ug_.begin(), get(graph::vertex_index, graph_));
		}


		auto A()
		{
			return make_iterator_property_map(A_.begin(), get(graph::edge_index, graph_));
		}


		auto A() const
		{
			return make_iterator_property_map(A_.begin(), get(graph::edge_index, graph_));
		}


		auto B()
		{
			return make_iterator_property_map(B_.begin(), get(graph::edge_index, graph_));
		}


		auto B() const
		{
			return make_iterator_property_map(B_.begin(), get(graph::edge_index, graph_));
		}


		auto b()
		{
			return make_iterator_property_map(b_.begin(), get(graph::edge_index, graph_));
		}


		auto b() const
		{
			return make_iterator_property_map(b_.begin(), get(graph::edge_index, graph_));
		}


		auto x() const
		{
			return make_iterator_property_map(x_.begin(), get(graph::vertex_index, graph_));
		}


		auto x()
		{
			return make_iterator_property_map(x_.begin(), get(graph::vertex_index, graph_));
		}


		auto u() const
		{
			return make_iterator_property_map(u_.begin(), get(graph::vertex_index, graph_));
		}


		auto u()
		{
			return make_iterator_property_map(u_.begin(), get(graph::vertex_index, graph_));
		}


		auto pi() const
		{
			return make_iterator_property_map(pi_.begin(), get(graph::edge_index, graph_));
		}

		
	private:
		OcpGraph graph_;
		std::vector<OcpSize> size_;

		// --------------------------------
		//
		// QP problem data
		//
		// --------------------------------
		std::vector<UnpaddedMatrix<Real, SO>> A_;
		std::vector<UnpaddedMatrix<Real, SO>> B_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> b_;
		std::vector<UnpaddedMatrix<Real, SO>> Q_;
		std::vector<UnpaddedMatrix<Real, SO>> S_;
		std::vector<UnpaddedMatrix<Real, SO>> R_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> q_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> r_;

		std::vector<UnpaddedMatrix<Real, SO>> Zl_;
		std::vector<UnpaddedMatrix<Real, SO>> Zu_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> zl_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> zu_;

		std::vector<UnpaddedMatrix<Real, SO>> C_;
		std::vector<UnpaddedMatrix<Real, SO>> D_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lg_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> ug_;

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
		
		// Lower and upper state and input bounds
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lx_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> ux_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lu_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> uu_;

		// Lambda multipliers for lower and upper state and input bounds
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_lx_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_ux_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_lu_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_uu_;

		// --------------------------------
		//
		// QP solution data
		//
		// --------------------------------

		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> x_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> u_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> ls_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> us_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> pi_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_lb_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_ub_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_lg_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_ug_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_ls_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lam_us_;

		// --------------------------------
		//
		// HPIPM solver data
		//
		// --------------------------------
		// OCP QP dimensions structure is intentionally stored on heap.
		// This is because we want the move ctor for HpipmWorkspace to work,
		// and othet structures that are members of HpipmWorkspace hold a pointer to OcpQpDim.
		std::unique_ptr<hpipm::OcpQpDim<Real>> ocpQpDim_;
		std::unique_ptr<hpipm::OcpQpIpmArg<Real>> solverArg_;
		std::unique_ptr<hpipm::OcpQp<Real>> qp_;
		std::unique_ptr<hpipm::OcpQpSol<Real>> sol_;
		std::unique_ptr<hpipm::OcpQpIpmWs<Real>> solverWorkspace_;


		// Warmstarting disabled on purpose.
		// On AMD K8 (hpmpc compiled for SSE3), WITHOUT warmstarting it is significantly
		// FASTER (9ms vs 14ms per time step) than with warmstarting. I am curious why.
		bool _warmStart = false;


		class InitVisitor
        :   public graph::default_bfs_visitor 
        {
        public:
            InitVisitor(HpipmWorkspace& ws)
            :   ws_{ws}
            {
            }

        
			template <typename Graph>
            void discover_vertex(OcpVertexDescriptor u, Graph const& g)
            {
				auto const& sz = get(ws_.size(), u);

				{
					// Init idxbx to 0,1,2,... s.t. maximum necessary workspace size is returned.
					std::unique_ptr<int []> idxbx {new int[sz.nx()]};
					for (int i = 0; i < sz.nx(); ++i)
						idxbx[i] = i;
					
					ws_.hidxbx_.emplace_back(std::move(idxbx));

					// Init idxbu to 0,1,2,... s.t. maximum necessary workspace size is returned.
					std::unique_ptr<int []> idxbu {new int [sz.nu()]};
					for (int i = 0; i < sz.nu(); ++i)
						idxbu[i] = i;

					ws_.hidxbu_.emplace_back(std::move(idxbu));
				}
					
				{
					// Initialize idxs to 0,1,2, ... .
					std::unique_ptr<int []> idxs {new int[sz.ns()]};
					for (int i = 0; i < sz.ns(); ++i)
						idxs[i] = i;

					ws_.idxs_.emplace_back(std::move(idxs));
				}

				ws_.Q_.emplace_back(sz.nx(), sz.nx());
				ws_.S_.emplace_back(sz.nu(), sz.nx());
				ws_.R_.emplace_back(sz.nu(), sz.nu());
				ws_.q_.emplace_back(sz.nx());
				ws_.r_.emplace_back(sz.nu());
				ws_.Zl_.emplace_back(sz.ns(), sz.ns());
				ws_.Zu_.emplace_back(sz.ns(), sz.ns());
				ws_.zl_.emplace_back(sz.ns());
				ws_.zu_.emplace_back(sz.ns());				

				ws_.lu_.emplace_back(sz.nu(), -inf<Real>());
				ws_.lx_.emplace_back(sz.nx(), -inf<Real>());

				ws_.uu_.emplace_back(sz.nu(), inf<Real>());
				ws_.ux_.emplace_back(sz.nx(), inf<Real>());
				
				ws_.lbls_.emplace_back(sz.ns());
				ws_.lbus_.emplace_back(sz.ns());

				ws_.C_ .emplace_back(sz.nc(), sz.nx());
				ws_.D_ .emplace_back(sz.nc(), sz.nu());
				ws_.lg_.emplace_back(sz.nc());
				ws_.ug_.emplace_back(sz.nc());

				ws_.x_.emplace_back(blaze::ZeroVector<Real>(sz.nx()));
				ws_.u_.emplace_back(blaze::ZeroVector<Real>(sz.nu()));
				ws_.ls_.emplace_back(sz.ns());
				ws_.us_.emplace_back(sz.ns());
				
				ws_.lam_lb_.emplace_back(sz.nx() + sz.nu());
				ws_.lam_ub_.emplace_back(sz.nx() + sz.nu());
				ws_.lam_lg_.emplace_back(sz.nc());
				ws_.lam_ug_.emplace_back(sz.nc());
				ws_.lam_ls_.emplace_back(sz.ns());
				ws_.lam_us_.emplace_back(sz.ns());

				ws_.lam_lx_.emplace_back(sz.nx(), sNaN<Real>());
				ws_.lam_ux_.emplace_back(sz.nx(), sNaN<Real>());
				ws_.lam_lu_.emplace_back(sz.nx(), sNaN<Real>());
				ws_.lam_uu_.emplace_back(sz.nx(), sNaN<Real>());
            }


			template <typename Graph>
			void tree_edge(OcpEdgeDescriptor e, Graph const& g)
			{
				auto const u = source(e, g);
				auto const v = target(e, g);

				if (v != u + 1)
					throw std::invalid_argument("Invalid tree structure in HpipmWorkspace ctor:	vertices are not sequentially connected");

				auto const& sz_u = get(ws_.size(), u);
				auto const& sz_v = get(ws_.size(), v);

				ws_.A_.emplace_back(sz_v.nx(), sz_u.nx());
				ws_.B_.emplace_back(sz_v.nx(), sz_u.nu());
				ws_.b_.emplace_back(sz_v.nx());
				ws_.pi_.emplace_back(sz_v.nx());
			}


			template <typename Graph>
			void non_tree_edge(OcpEdgeDescriptor e, Graph const& g)
			{
				throw std::invalid_argument("Invalid tree structure in HpipmWorkspace ctor:	non-tree graph structure detected.");
			}

        
        private:
            HpipmWorkspace& ws_;
        };


		template <typename SizeMap>
		static std::unique_ptr<hpipm::OcpQpDim<Real>> makeOcpQpDim(OcpGraph const& graph, SizeMap size)
		{
			auto const n = num_vertices(graph);
			if (n == 0)
				TMPC_THROW_EXCEPTION(std::invalid_argument("Graph is empty"));

			auto dim = std::make_unique<hpipm::OcpQpDim<Real>>(n - 1);
			for (auto v : graph::vertices(graph))
			{
				dim->set_nx(v, get(size, v).nx());
				dim->set_nu(v, get(size, v).nu());
				dim->set_nbx(v, get(size, v).nx());
				dim->set_nbu(v, get(size, v).nu());
				dim->set_ng(v, get(size, v).nc());
				dim->set_nsbx(v, 0);
				dim->set_nsbu(v, 0);
				dim->set_nsg(v, 0);
			}

			return dim;
		}


		/// @brief Re-create solver workspace with new args
		void recreateSolverWorkspace()
		{
			solverWorkspace_.reset(new hpipm::OcpQpIpmWs<Real> {*ocpQpDim_, *solverArg_});
		}
	};


	template <typename Real>
    struct RealOf<HpipmWorkspace<Real>>
    {
        using type = Real;
    };
}
