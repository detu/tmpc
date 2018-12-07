#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/core/PropertyMap.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_d_ocp_qp_kkt.h>

#include "detail/VectorPtrPropertyMap.hpp"
#include "detail/MatrixPtrPropertyMap.hpp"

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>


namespace tmpc
{
	namespace detail
	{
		template <typename Real>
		struct HpipmApi;

		template <>
		struct HpipmApi<double>
		{
			using ocp_qp = ::d_ocp_qp;
			using ocp_qp_sol = ::d_ocp_qp_sol;
			using ocp_qp_ipm_workspace = ::d_ocp_qp_ipm_workspace;
			using ocp_qp_ipm_arg = ::d_ocp_qp_ipm_arg;
			using ocp_qp_dim = ::d_ocp_qp_dim;

			static auto constexpr memsize_ocp_qp = ::d_memsize_ocp_qp;
			static auto constexpr create_ocp_qp = ::d_create_ocp_qp;
			static auto constexpr cvt_colmaj_to_ocp_qp = ::d_cvt_colmaj_to_ocp_qp;
			static auto constexpr memsize_ocp_qp_sol = ::d_memsize_ocp_qp_sol;
			static auto constexpr create_ocp_qp_sol = ::d_create_ocp_qp_sol;
			static auto constexpr cvt_ocp_qp_sol_to_colmaj = ::d_cvt_ocp_qp_sol_to_colmaj;
			static auto constexpr memsize_ocp_qp_ipm = ::d_memsize_ocp_qp_ipm;
			static auto constexpr create_ocp_qp_ipm = ::d_create_ocp_qp_ipm;
			static auto constexpr solve_ocp_qp_ipm = ::d_solve_ocp_qp_ipm;
		};
	}


	class HpipmException
	:	public std::runtime_error
	{
	public:
		HpipmException(int code);
		int code() const { return _code; }

	private:
		int const _code;
	};


	/**
	 * \brief Multistage QP solver using HPIPM
	 *
	 * \tparam <Kernel_> matrix kernel type
	 */
	template <typename Kernel_>
	class HpipmWorkspace
	{
	public:
		using Kernel = Kernel_;
		static StorageOrder constexpr SO = StorageOrder::columnMajor;
		using Real = typename Kernel::Real;
		

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
		,	ocpQpDim_ {}
		,	solverArg_ {}
		{
			copyProperty(size, iterator_property_map(size_.begin(), get(graph::vertex_index, graph_)), 
				graph::vertices(graph_));
			
			auto const nv = num_vertices(graph_);
			auto const ne = num_edges(graph_);

			if (nv < 2)
				throw std::invalid_argument("HPIPM needs an at least 2-stages problem");

			if (ne + 1 != num_vertices(graph_))
				throw std::invalid_argument("Number of edges in HPIPM size graph must be 1 less than the number of vertices");

			solverArg_.alpha_min = 1e-8;
			solverArg_.res_g_max = 1e-8;
			solverArg_.res_b_max = 1e-8;
			solverArg_.res_d_max = 1e-12;
			solverArg_.res_m_max = 1e-12;
			solverArg_.mu0 = 2.0;
			solverArg_.iter_max = max_iter;
			solverArg_.stat_max = 100;
			solverArg_.pred_corr = 1;
			
			// Calculate total number of int and Real elements to be allocated.
			size_t count_real = 0, count_int = 0;
			breadth_first_search(graph_, vertex(0, graph_), visitor(ElementCountVisitor(this->size(), count_real, count_int)));

			// Allocate memory pools
			realPool_.resize(count_real, sNaN<Real>());
			intPool_.resize(count_int, -1);
			
			// Preallocate arrays holding QP vertex data pointers.
			nx_.reserve(nv);
			nu_.reserve(nv);
			nb_.reserve(nv);
			nbx_.reserve(nv);
			nbu_.reserve(nv);
			ng_.reserve(nv);
			ns_.reserve(nv);
			nsbx_.reserve(nv);
			nsbu_.reserve(nv);
			nsg_.reserve(nv);
			lbls_.reserve(nv);
			lbus_.reserve(nv);
			hidxb_.reserve(nv);
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
			lb_.reserve(nv);
			ub_.reserve(nv);
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

			ocpQpDim_.N = nv - 1;
			ocpQpDim_.nx = nx_.data();
			ocpQpDim_.nu = nu_.data();
			ocpQpDim_.nb = nb_.data();
			ocpQpDim_.nbx = nbx_.data();
			ocpQpDim_.nbu = nbu_.data();
			ocpQpDim_.ng = ng_.data();
			ocpQpDim_.ns = ns_.data();
			ocpQpDim_.nsbx = nsbx_.data();
			ocpQpDim_.nsbu = nsbu_.data();
			ocpQpDim_.nsg = nsg_.data();
				
			if (nv > 1)
			{
				// Allocate big enough memory pools for Qp, QpSol and SolverWorkspace.
				// nb_ is set to nx+nu at this point, which ensures maximum capacity.
				ocpQpMem_.resize(Hpipm::memsize_ocp_qp(&ocpQpDim_));
				ocpQpSolMem_.resize(Hpipm::memsize_ocp_qp_sol(&ocpQpDim_));
				resizeSolverWorkspace();			
			}
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
	
				// Call HPIPM
				auto const ret = Hpipm::solve_ocp_qp_ipm(&qp, &sol, &solverArg_, &solver_workspace);
				numIter_ = solver_workspace.iter;
	
				if (ret != 0)
					throw HpipmException(ret);

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
			d_fact_solve_kkt_unconstr_ocp_qp(&qp, &sol, &solverArg_, &solver_workspace);
		}


		/// @brief Copy QP to HPIPM data structures.
		void convertProblem()
		{
			updateBounds();
			
			// Init HPIPM structures. Since the nb_, nbx_, nbu_, might change if the bounds change, 
			// we need to do it every time, sorry.
			// Hopefully these operations are not expensive compared to actual solving.
			Hpipm::create_ocp_qp(&ocpQpDim_, &qp, ocpQpMem_.data());
			Hpipm::create_ocp_qp_sol(&ocpQpDim_, &sol, ocpQpSolMem_.data());
			Hpipm::create_ocp_qp_ipm(&ocpQpDim_, &solverArg_, &solver_workspace, solverWorkspaceMem_.data());

			// Convert the problem
			Hpipm::cvt_colmaj_to_ocp_qp(
				A_.data(), B_.data(), b_.data(), 
				Q_.data(), S_.data(), R_.data(), q_.data(), r_.data(), 
				hidxb_.data(), lb_.data(), ub_.data(), 
				C_.data(), D_.data(), lg_.data(), ug_.data(), 
				Zl_.data(), Zu_.data(), zl_.data(), zu_.data(), idxs_.data(),
				lbls_.data(), lbus_.data(),
				&qp);
		}


		/// @brief Copy QP solution from HPIPM data structures.
		void convertSolution()
		{
			// Convert the solution
			Hpipm::cvt_ocp_qp_sol_to_colmaj(&sol, 
				u_.data(), x_.data(), ls_.data(), us_.data(),
				pi_.data(), lam_lb_.data(), lam_ub_.data(), lam_lg_.data(), lam_ug_.data(),
				lam_ls_.data(), lam_us_.data());
		}


		size_t maxIter() const noexcept 
		{ 
			return solverArg_.iter_max; 
		}

		
		void maxIter(size_t val) 
		{ 
			if (val != solverArg_.iter_max)
			{
				solverArg_.iter_max = val;
				resizeSolverWorkspace();
			}
		}


		/// \brief Get number of iterations performed by the QP solver.
		unsigned numIter() const { return numIter_; }


		void alphaMin(Real val)
		{
			if (val <= 0.)
				throw std::invalid_argument("HpipmWorkspace::alphaMin(): value must be positive");

			solverArg_.alpha_min = val;
		}

		
		void resGMax(Real val)
		{
			if (val <= 0.)
				throw std::invalid_argument("HpipmWorkspace::resGMax(): value must be positive");

			solverArg_.res_g_max = val;
		}

		
		void resBMax(Real val)
		{
			if (val <= 0.)
				throw std::invalid_argument("HpipmWorkspace::resBMax(): value must be positive");
			
			solverArg_.res_b_max = val;
		}

		
		void resDMax(Real val)
		{
			if (val <= 0.)
				throw std::invalid_argument("HpipmWorkspace::resDMax(): value must be positive");
			
			solverArg_.res_d_max = val;
		}

		
		void resMMax(Real val)
		{
			if (val <= 0.)
				throw std::invalid_argument("HpipmWorkspace::resMMax(): value must be positive");
		
			solverArg_.res_m_max = val;
		}


		void mu0(Real val)
		{
			if (val <= 0.)
				throw std::invalid_argument("HpipmWorkspace::mu0(): value must be positive");
		
			solverArg_.mu0 = val;
		}
		
		
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
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(Q_.begin(), get(graph::vertex_index, graph_)), size_Q(size()));
		}


		auto Q() const
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(Q_.begin(), get(graph::vertex_index, graph_)), size_Q(size()));
		}


		auto R()
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(R_.begin(), get(graph::vertex_index, graph_)), size_R(size()));
		}


		auto R() const
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(R_.begin(), get(graph::vertex_index, graph_)), size_R(size()));
		}


		auto S()
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(S_.begin(), get(graph::vertex_index, graph_)), size_S(size()));
		}


		auto S() const
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(S_.begin(), get(graph::vertex_index, graph_)), size_S(size()));
		}


		auto q()
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(q_.begin(), get(graph::vertex_index, graph_)), size_x(size()));
		}


		auto q() const
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(q_.begin(), get(graph::vertex_index, graph_)), size_x(size()));
		}


		auto r()
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(r_.begin(), get(graph::vertex_index, graph_)), size_u(size()));
		}


		auto r() const
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(r_.begin(), get(graph::vertex_index, graph_)), size_u(size()));
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
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(C_.begin(), get(graph::vertex_index, graph_)), size_C(size()));
		}


		auto C() const
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(C_.begin(), get(graph::vertex_index, graph_)), size_C(size()));
		}


		auto D()
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(D_.begin(), get(graph::vertex_index, graph_)), size_D(size()));
		}


		auto D() const
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(D_.begin(), get(graph::vertex_index, graph_)), size_D(size()));
		}


		auto ld()
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(lg_.begin(), get(graph::vertex_index, graph_)), size_d(size()));
		}


		auto ld() const
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(lg_.begin(), get(graph::vertex_index, graph_)), size_d(size()));
		}


		auto ud()
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(ug_.begin(), get(graph::vertex_index, graph_)), size_d(size()));
		}


		auto ud() const
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(ug_.begin(), get(graph::vertex_index, graph_)), size_d(size()));
		}


		auto A()
		{
			return detail::makeMatrixPtrPropertyMap<OcpEdgeDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(A_.begin(), get(graph::edge_index, graph_)), size_A(size(), graph_));
		}


		auto A() const
		{
			return detail::makeMatrixPtrPropertyMap<OcpEdgeDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(A_.begin(), get(graph::edge_index, graph_)), size_A(size(), graph_));
		}


		auto B()
		{
			return detail::makeMatrixPtrPropertyMap<OcpEdgeDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(B_.begin(), get(graph::edge_index, graph_)), size_B(size(), graph_));
		}


		auto B() const
		{
			return detail::makeMatrixPtrPropertyMap<OcpEdgeDescriptor, CustomMatrix<Kernel, unaligned, unpadded, SO>>(
				make_iterator_property_map(B_.begin(), get(graph::edge_index, graph_)), size_B(size(), graph_));
		}


		auto b()
		{
			return detail::makeVectorPtrPropertyMap<OcpEdgeDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(b_.begin(), get(graph::edge_index, graph_)), size_b(size(), graph_));
		}


		auto b() const
		{
			return detail::makeVectorPtrPropertyMap<OcpEdgeDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(b_.begin(), get(graph::edge_index, graph_)), size_b(size(), graph_));
		}


		auto x() const
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(x_.begin(), get(graph::vertex_index, graph_)), size_x(size()));
		}


		auto u() const
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(u_.begin(), get(graph::vertex_index, graph_)), size_u(size()));
		}


		auto pi() const
		{
			return detail::makeVectorPtrPropertyMap<OcpEdgeDescriptor, CustomVector<Kernel, unaligned, unpadded>>(
				make_iterator_property_map(pi_.begin(), get(graph::edge_index, graph_)), size_b(size(), graph_));
		}

		
	private:
		using Hpipm = detail::HpipmApi<Real>;
		
		OcpGraph graph_;
		std::vector<OcpSize> size_;

		std::vector<Real> realPool_;
		std::vector<int> intPool_;

		// --------------------------------
		//
		// QP problem data
		//
		// --------------------------------
		std::vector<Real *> A_;
		std::vector<Real *> B_;
		std::vector<Real *> b_;
		std::vector<Real *> Q_;
		std::vector<Real *> S_;
		std::vector<Real *> R_;
		std::vector<Real *> q_;
		std::vector<Real *> r_;

		std::vector<Real *> Zl_;
		std::vector<Real *> Zu_;
		std::vector<Real *> zl_;
		std::vector<Real *> zu_;

		std::vector<Real *> lb_;
		std::vector<Real *> ub_;
		std::vector<Real *> C_;
		std::vector<Real *> D_;
		std::vector<Real *> lg_;
		std::vector<Real *> ug_;

		// Array of NX sizes
		std::vector<int> nx_;

		// Array of NU sizes
		std::vector<int> nu_;

		// Array of NB (bound constraints) sizes
		std::vector<int> nb_;

		// Array of NBX (state bound constraints) sizes
		std::vector<int> nbx_;

		// Array of NBU (input bound constraints) sizes
		std::vector<int> nbu_;

		// Array of NG (path constraints) sizes
		std::vector<int> ng_;

		// Array of NS (soft constraints) sizes
		std::vector<int> ns_;

		// Array of NSBX (soft state bound constraints) sizes
		std::vector<int> nsbx_;

		// Array of NSBU (soft input bound constraints) sizes
		std::vector<int> nsbu_;

		// Array of NSG (soft general bound constraints) sizes
		std::vector<int> nsg_;

		// Lower bound of lower slack.
		std::vector<Real *> lbls_;

		// Upper bound of upper slack.
		std::vector<Real *> lbus_;

		// Hard constraints index
		std::vector<int *> hidxb_;

		// Soft constraints index
		std::vector<int *> idxs_;
		
		// Full lower and upper state and input bounds (can contain +-infs)
		std::vector<DynamicVector<Kernel>> lx_;
		std::vector<DynamicVector<Kernel>> ux_;
		std::vector<DynamicVector<Kernel>> lu_;
		std::vector<DynamicVector<Kernel>> uu_;

		// Full lambda multipliers for lower and upper state and input bounds
		std::vector<DynamicVector<Kernel>> lam_lx_;
		std::vector<DynamicVector<Kernel>> lam_ux_;
		std::vector<DynamicVector<Kernel>> lam_lu_;
		std::vector<DynamicVector<Kernel>> lam_uu_;

		// --------------------------------
		//
		// QP solution data
		//
		// --------------------------------

		std::vector<Real *> x_;
		std::vector<Real *> u_;
		std::vector<Real *> ls_;
		std::vector<Real *> us_;
		std::vector<Real *> pi_;
		std::vector<Real *> lam_lb_;
		std::vector<Real *> lam_ub_;
		std::vector<Real *> lam_lg_;
		std::vector<Real *> lam_ug_;
		std::vector<Real *> lam_ls_;
		std::vector<Real *> lam_us_;

		/// \brief Number of iterations performed by the QP solver.
		int numIter_ = 0;

		// --------------------------------
		//
		// HPIPM solver data
		//
		// --------------------------------
		std::vector<char> ocpQpMem_;
		std::vector<char> ocpQpSolMem_;
		std::vector<char> solverWorkspaceMem_;

		typename Hpipm::ocp_qp_dim ocpQpDim_;
		typename Hpipm::ocp_qp_ipm_arg solverArg_;
		typename Hpipm::ocp_qp qp {};
		typename Hpipm::ocp_qp_sol sol {};
		typename Hpipm::ocp_qp_ipm_workspace solver_workspace {};


		// Warmstarting disabled on purpose.
		// On AMD K8 (hpmpc compiled for SSE3), WITHOUT warmstarting it is significantly
		// FASTER (9ms vs 14ms per time step) than with warmstarting. I am curious why.
		bool _warmStart = false;


		// This visitor calculates the total number of elements needed
		// for all arrays passed to and returned from HPMPC.
		template <typename SizeMap>
		class ElementCountVisitor
        :   public graph::default_bfs_visitor 
        {
        public:
            ElementCountVisitor(SizeMap size_map, size_t& real_count, size_t& int_count)
            :   sizeMap_{size_map}
			,	realCount_{real_count}
			,	intCount_{int_count}
            {
				realCount_ = 0;
				intCount_ = 0;
            }

        
			template <typename Graph>
            void discover_vertex(OcpVertexDescriptor u, Graph const& g)
            {
				auto const& sz = get(sizeMap_, u);
				auto const max_nb = sz.nu() + sz.nx();

				intCount_ += max_nb;	// idxb
				intCount_ += sz.ns();	// idxs
                
				realCount_ += sz.nx() * sz.nx();	// Q
				realCount_ += sz.nu() * sz.nx();	// S
				realCount_ += sz.nu() * sz.nu();	// R 
				realCount_ += sz.nx();	// q
				realCount_ += sz.nu();	// r
				realCount_ += sz.ns() * sz.ns();	// Zl
				realCount_ += sz.ns() * sz.ns();	// Zu
				realCount_ += sz.ns();	// zl
				realCount_ += sz.ns();	// zu
				realCount_ += max_nb;	// lb
				realCount_ += max_nb;	// ub
				realCount_ += sz.ns();	// lbls
				realCount_ += sz.ns();	// lbus

				realCount_ += sz.nc() * sz.nx();	// C
				realCount_ += sz.nc() * sz.nu();	// D
				realCount_ += sz.nc();	// lg
				realCount_ += sz.nc();	// ug

				realCount_ += sz.nx();	// x
				realCount_ += sz.nu();	// u
				realCount_ += sz.ns();	// ls
				realCount_ += sz.ns();	// us
				
				realCount_ += max_nb;	// lam_lb
				realCount_ += max_nb;	// lam_ub
				realCount_ += sz.nc();	// lam_lg
				realCount_ += sz.nc();	// lam_ug
				
				realCount_ += sz.ns();	// lam_ls
				realCount_ += sz.ns();	// lam_us
            }


			template <typename Graph>
			void tree_edge(OcpEdgeDescriptor e, Graph const& g)
			{
				auto const& sz_u = get(sizeMap_, source(e, g));
				auto const& sz_v = get(sizeMap_, target(e, g));

				realCount_ += sz_v.nx() * sz_u.nx();	// A
				realCount_ += sz_v.nx() * sz_u.nu();	// B
				realCount_ += sz_v.nx();	// b
				realCount_ += sz_v.nx();	// pi
			}

        
        private:
            SizeMap sizeMap_;
			size_t& realCount_;
			size_t& intCount_;
        };


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
				auto const max_nb = sz.nu() + sz.nx();

				ws_.nx_.push_back(sz.nx());
				ws_.nu_.push_back(sz.nu());
				ws_.nb_.push_back(max_nb);	// upper estimate
				ws_.nbx_.push_back(sz.nx());	// upper estimate
				ws_.nbu_.push_back(sz.nu());	// upper estimate
				ws_.ng_.push_back(sz.nc());
				ws_.ns_.push_back(sz.ns());
				ws_.nsbx_.push_back(sz.nx());	// upper estimate
				ws_.nsbu_.push_back(sz.nu());	// upper estimate
				ws_.nsg_.push_back(sz.nc());	// upper estimate
				
				{
					// Allocate idxb array. It will be properly initialized in solve().
					int * const idxb = allocInt(max_nb);
					ws_.hidxb_.push_back(idxb);

					// Init idxb to 0,1,2,... s.t. maximum necessary workspace size is returned.
					for (int i = 0; i < max_nb; ++i)
						idxb[i] = i;
				}
					
				{
					// Allocate idxs array.				
					int * const idxs = allocInt(sz.ns());
					ws_.idxs_.push_back(idxs);
				
					// Initialize idxs to 0,1,2, ... .
					for (int i = 0; i < sz.ns(); ++i)
						idxs[i] = i;
				}

				ws_.Q_.push_back(allocReal(sz.nx() * sz.nx()));
				ws_.S_.push_back(allocReal(sz.nu() * sz.nx()));
				ws_.R_.push_back(allocReal(sz.nu() * sz.nu()));
				ws_.q_.push_back(allocReal(sz.nx()));
				ws_.r_.push_back(allocReal(sz.nu()));
				ws_.Zl_.push_back(allocReal(sz.ns() * sz.ns()));
				ws_.Zu_.push_back(allocReal(sz.ns() * sz.ns()));
				ws_.zl_.push_back(allocReal(sz.ns()));
				ws_.zu_.push_back(allocReal(sz.ns()));				

				ws_.lb_.push_back(allocReal(max_nb));	// the array will be updated during solve()
				ws_.lu_.emplace_back(sz.nu(), -inf<Real>());
				ws_.lx_.emplace_back(sz.nx(), -inf<Real>());

				ws_.ub_.push_back(allocReal(max_nb));	// the array will be updated during solve()
				ws_.uu_.emplace_back(sz.nu(), inf<Real>());			
				ws_.ux_.emplace_back(sz.nx(), inf<Real>());
				
				ws_.lbls_.push_back(allocReal(sz.ns()));
				ws_.lbus_.push_back(allocReal(sz.ns()));

				ws_.C_ .push_back(allocReal(sz.nc() * sz.nx()));
				ws_.D_ .push_back(allocReal(sz.nc() * sz.nu()));
				ws_.lg_.push_back(allocReal(sz.nc()));
				ws_.ug_.push_back(allocReal(sz.nc()));

				ws_.x_.push_back(allocReal(sz.nx()));
				ws_.u_.push_back(allocReal(sz.nu()));
				ws_.ls_.push_back(allocReal(sz.ns()));
				ws_.us_.push_back(allocReal(sz.ns()));
				
				ws_.lam_lb_.push_back(allocReal(max_nb));
				ws_.lam_ub_.push_back(allocReal(max_nb));
				ws_.lam_lg_.push_back(allocReal(sz.nc()));
				ws_.lam_ug_.push_back(allocReal(sz.nc()));
				ws_.lam_ls_.push_back(allocReal(sz.ns()));
				ws_.lam_us_.push_back(allocReal(sz.ns()));

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

				ws_.A_.push_back(allocReal(sz_v.nx() * sz_u.nx()));
				ws_.B_.push_back(allocReal(sz_v.nx() * sz_u.nu()));
				ws_.b_.push_back(allocReal(sz_v.nx()));
				ws_.pi_.push_back(allocReal(sz_v.nx()));
			}


			template <typename Graph>
			void non_tree_edge(OcpEdgeDescriptor e, Graph const& g)
			{
				throw std::invalid_argument("Invalid tree structure in HpipmWorkspace ctor:	non-tree graph structure detected.");
			}

        
        private:
            HpipmWorkspace& ws_;
			size_t allocatedReal_ = 0;
			size_t allocatedInt_ = 0;


			Real * allocReal(size_t n)
			{
				if (allocatedReal_ + n > ws_.realPool_.size())
					throw std::bad_alloc();

				Real * const ptr = ws_.realPool_.data() + allocatedReal_;
				allocatedReal_ += n;

				return ptr;
			}


			int * allocInt(size_t n)
			{
				if (allocatedInt_ + n > ws_.intPool_.size())
					throw std::bad_alloc();

				int * const ptr = ws_.intPool_.data() + allocatedInt_;
				allocatedInt_ += n;

				return ptr;
			}
        };


		class PopulateBounds
		{
		public:
			PopulateBounds(Real * lb, Real * ub, int * idxb)
			:	lb_(lb)
			,	ub_(ub)
			,	idxb_(idxb)
			{
			}


			bool operator()(Real l, Real u)
			{
				bool bounded = false;
				
				if (std::isfinite(l) && std::isfinite(u))
				{
					// If both bounds are finite, add i to the bounds index,
					// and copy values to the lb_internal_ and ub_internal_.
					*idxb_++ = count_;
					*lb_++ = l;
					*ub_++ = u;
					bounded = true;
				}
				else 
				{
					// Otherwise, check that the values are [-inf, inf]
					if (!(l == -inf<Real>() && u == inf<Real>()))
						throw std::invalid_argument("An invalid QP bound is found. For HPMPC/HPIPM, "
							"the bounds should be either both finite or [-inf, inf]");
				}

				++count_;
				
				return bounded;
			}


		private:
			Real * lb_;
			Real * ub_;
			int * idxb_;
			int count_ = 0;
		};


		void resizeSolverWorkspace()
		{
			solverWorkspaceMem_.resize(Hpipm::memsize_ocp_qp_ipm(&ocpQpDim_, &solverArg_));
		}


		/// @brief Recalculate bounds indices of each stage, as the bounds might have changed.
		void updateBounds()
		{
			auto const vertex_id = get(graph::vertex_index, graph_);
				
			// Update the nb_ array.
			for(auto v : graph::vertices(graph_))
			{
				auto const v_i = get(vertex_id, v);
				auto const& sz = get(size(), v);

				int nbx = 0, nbu = 0;
				PopulateBounds pb(lb_[v_i], ub_[v_i], hidxb_[v_i]);
				
				// Cycle through the bounds and check for infinities
				{
					auto const& lu = lu_[v_i];
					auto const& uu = uu_[v_i];
					
					for (size_t i = 0; i < sz.nu(); ++i)
						if (pb(lu[i], uu[i]))
							++nbu;
				}

				{
					auto const& lx = lx_[v_i];
					auto const& ux = ux_[v_i];
				
					for (size_t i = 0; i < sz.nx(); ++i)
						if (pb(lx[i], ux[i]))
							++nbx;
				}
					
				nb_[v_i] = nbx + nbu;
				nbx_[v_i] = nbx;
				nbu_[v_i] = nbu;
			
				// Cycle through the soft bounds and count number of soft state bounds (nsbx),
				// number of soft input bounds (nsbu) and number of soft general bounds (nsg).
				int nsbx = 0, nsbu = 0, nsg = 0;
				int const * const idxs = idxs_[v_i];
				
				for (size_t i = 0; i < sz.ns(); ++i)
				{
					if (sz.nu() <= idxs[i] && idxs[i] < sz.nu() + sz.nx())
						++nsbx;
						
					if (idxs[i] < sz.nu())
						++nsbu;
						
					if (idxs[i] >= sz.nu() + sz.nx())
						++nsg;
				}
				
				nsbx_[v_i] = nsbx;
				nsbu_[v_i] = nsbu;
				nsg_[v_i] = nsg;
			}
		}
	};
}
