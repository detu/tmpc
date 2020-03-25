#pragma once

#include <tmpc/qp/QpSolverException.hpp>
#include <tmpc/qp/HpmpcErrorInfo.hpp>
//#include <tmpc/property_map/BundlePropertyMap.hpp>

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSolutionBase.hpp>
#include <tmpc/qp/OcpQpBase.hpp>
//#include <tmpc/qp/QpWorkspaceBase.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

//#include "detail/HpxxxVertexData.hpp"
//#include "detail/HpxxxEdgeData.hpp"
#include <tmpc/property_map/VectorPtrPropertyMap.hpp>
#include <tmpc/property_map/MatrixPtrPropertyMap.hpp>
//#include "TreeQpWorkspaceAdaptor.hpp"
//#include "OcpQpVertexElement.hpp"
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/Exception.hpp>
#include <tmpc/Traits.hpp>

#include <c_interface.h>

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <algorithm>
#include <cmath>
#include <vector>
#include <array>
#include <memory>

namespace tmpc
{
	// Type-parameterized HPMPC interface

	namespace hpmpc
	{
		inline int c_order_ip_ocp_hard_tv(
			int *kk, int k_max, double mu0, double mu_tol,
			int N, int const *nx, int const *nu, int const *nb, int const * const *hidxb, int const *ng, int N2,
			int warm_start,
			double const * const *A, double const * const *B, double const * const *b,
			double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r,
			double const * const *lb, double const * const *ub,
			double const * const *C, double const * const *D, double const * const *lg, double const * const *ug,
			double * const *x, double * const *u, double * const *pi, double * const *lam,
			double *inf_norm_res,
			void *work0,
			double *stat)
		{
			return ::c_order_d_ip_ocp_hard_tv(
				kk, k_max, mu0, mu_tol,
				N, const_cast<int*>(nx), const_cast<int*>(nu), const_cast<int*>(nb), const_cast<int **>(hidxb), const_cast<int*>(ng), N2,
				warm_start,
				const_cast<double**>(A), const_cast<double**>(B), const_cast<double**>(b),
				const_cast<double**>(Q), const_cast<double**>(S), const_cast<double**>(R), const_cast<double**>(q), const_cast<double**>(r),
				const_cast<double**>(lb), const_cast<double**>(ub),
				const_cast<double**>(C), const_cast<double**>(D), const_cast<double**>(lg), const_cast<double**>(ug),
				const_cast<double**>(x), const_cast<double**>(u), const_cast<double**>(pi), const_cast<double**>(lam), 
				inf_norm_res,
				work0,
				stat);
		}

		
		template <typename Real>
		int ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int const * const * hidxb, int const *ng, int N2);


		template <>
		inline int ip_ocp_hard_tv_work_space_size_bytes<double>(int N, int const *nx, int const *nu, int const *nb, int const * const * hidxb, int const *ng, int N2)
		{
			return ::hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(
				N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), const_cast<int **>(hidxb), const_cast<int*>(ng), N2);
		}
	}


	/**
	 * \brief Multistage QP solver using HPMPC
	 *
	 * \tparam <Kernel> the math kernel type
	 */
	template <typename Real>
	class HpmpcWorkspace
	{
	public:
		static StorageOrder constexpr SO = StorageOrder::rowMajor;
		
		// Iteration statistics. HPMPC returns 5 Real numbers per iteration:
		// step length for predictor and corrector, 
		// centering parameter, duality measure for predictor and corrector.
		using IterStat = std::array<Real, 5>;


		std::string impl_solverName() const
		{
			return "HPMPC";
		}

	
		auto stat() const
		{
			return boost::make_iterator_range(stat_.begin(), stat_.begin() + numIter_);
		}


		/**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <typename SizeMap>
		explicit HpmpcWorkspace(OcpGraph const& graph, SizeMap size, unsigned max_iter = 100)
		:	graph_ {graph}
		,	size_(num_vertices(graph))
		,	stat_ {max_iter}
		{
			copyProperty(size, iterator_property_map(size_.begin(), get(graph::vertex_index, graph_)), graph::vertices(graph_));
			std::fill(infNormRes_.begin(), infNormRes_.end(), sNaN<Real>());

			auto const nv = num_vertices(graph_);
			auto const ne = num_edges(graph_);

			if (nv < 2)
				TMPC_THROW_EXCEPTION(std::invalid_argument("HPMPC needs an at least 2-stages problem"));

			if (ne + 1 != num_vertices(graph_))
				TMPC_THROW_EXCEPTION(std::invalid_argument("Number of edges in HPMPC size graph must be 1 less than the number of vertices"));

			// Calculate total number of int and Real elements to be allocated.
			size_t count_real = 0, count_int = 0;
			breadth_first_search(graph_, vertex(0, graph_), visitor(ElementCountVisitor<decltype(this->size())>(this->size(), count_real, count_int)));

			// Allocate memory pools
			realPool_.resize(count_real, sNaN<Real>());
			intPool_.resize(count_int, -1);

			// Preallocate arrays holding QP vertex data pointers.
			nx_.reserve(nv);
			nu_.reserve(nv);
			nb_.reserve(nv);
			ng_.reserve(nv);
			hidxb_.reserve(nv);
	
			Q_ .reserve(nv);
			S_ .reserve(nv);
			R_ .reserve(nv);
			q_ .reserve(nv);
			r_ .reserve(nv);
			lb_.reserve(nv);
			ub_.reserve(nv);
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
			lam_.reserve(nv);

			// Preallocate arrays holding QP edge data pointers.
			A_ .reserve(ne);
			B_ .reserve(ne);
			b_ .reserve(ne);
			pi_.reserve(ne);

			// Traverse the graph and init vertex and edge data pointers.
			breadth_first_search(graph_, vertex(0, graph_), visitor(InitVisitor(*this)));

			// Allocate HPMPC working memory
			allocateSolverWorkspace();
		}


		/**
		 * \brief Copy constructor
		 *
		 * Copying is not allowed.
		 */
		HpmpcWorkspace(HpmpcWorkspace const&) = delete;


		/**
		 * \brief Move constructor
		 *
		 * Move-construction is ok.
		 */
		HpmpcWorkspace(HpmpcWorkspace&& rhs) = default;


		HpmpcWorkspace& operator= (HpmpcWorkspace const&) = delete;
		HpmpcWorkspace& operator= (HpmpcWorkspace &&) = delete;


		void solve()
		{
			// Recalculate bounds indices of each vertex, as the bounds might have changed.				
			// Update lb_, ub_, hidxb_, nb_.
			for(auto v : graph::vertices(graph_))
			{
				auto const v_i = get(graph::vertex_index, graph_, v);
				auto const& sz = get(size(), v);

				auto const& lx = lx_[v_i];
				auto const& ux = ux_[v_i];
				auto const& lu = lu_[v_i];
				auto const& uu = uu_[v_i];

				PopulateBounds pb(lb_[v_i], ub_[v_i], hidxb_[v_i], nb_[v_i]);

				// Cycle through the bounds and check for infinities
				for (size_t i = 0; i < sz.nu(); ++i)
					pb(lu[i], uu[i]);

				for (size_t i = 0; i < sz.nx(); ++i)
					pb(lx[i], ux[i]);
			}

			// Number of QP steps for HPMPC
			auto const N = num_vertices(graph_) - 1;

			// What is a good value for mu0?
			Real mu0 = 1.;

			// Call HPMPC
			auto const ret = hpmpc::c_order_ip_ocp_hard_tv(&numIter_, maxIter(), mu0, muTol_, N,
					nx_.data(), nu_.data(), nb_.data(), hidxb_.data(), ng_.data(), N, _warmStart ? 1 : 0, A_.data(), B_.data(), b_.data(),
					Q_.data(), S_.data(), R_.data(), q_.data(), r_.data(), lb_.data(), ub_.data(), C_.data(), D_.data(),
					lg_.data(), ug_.data(), x_.data(), u_.data(), pi_.data(), lam_.data(), infNormRes_.data(),
					solverWorkspace_.data(), stat_[0].data());

			if (ret != 0)
				TMPC_THROW_EXCEPTION(QpSolverException {} << HpmpcErrorInfo {ret});
		}


		std::size_t maxIter() const noexcept 
		{ 
			return stat_.size(); 
		}


		/// \brief Set max number of iterations
		void maxIter(size_t val)
		{
			stat_.resize(val);
		}


		Real muTol() const noexcept 
		{ 
			return muTol_; 
		}


		void muTol(Real val)
		{
			if (val <= 0.)
				throw std::invalid_argument("mu tolerance for hpmpc must be positive");

			muTol_ = val;
		}


		/// \brief Get number of iterations performed by the QP solver.
		auto numIter() const noexcept
		{ 
			return numIter_; 
		}


		auto size() const
		{
			return iterator_property_map(size_.begin(), get(graph::vertex_index, graph_));
		}


		auto const& graph() const noexcept
		{
			return graph_;
		}


		auto Q()
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(Q_.begin(), get(graph::vertex_index, graph_)), size_Q(size()));
		}


		auto Q() const noexcept
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(Q_.begin(), get(graph::vertex_index, graph_)), size_Q(size()));
		}


		auto R()
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(R_.begin(), get(graph::vertex_index, graph_)), size_R(size()));
		}


		auto R() const noexcept
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(R_.begin(), get(graph::vertex_index, graph_)), size_R(size()));
		}


		auto S()
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(S_.begin(), get(graph::vertex_index, graph_)), size_S(size()));
		}


		auto S() const noexcept
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(S_.begin(), get(graph::vertex_index, graph_)), size_S(size()));
		}


		auto C()
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(C_.begin(), get(graph::vertex_index, graph_)), size_C(size()));
		}


		auto C() const noexcept
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(C_.begin(), get(graph::vertex_index, graph_)), size_C(size()));
		}


		auto D()
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(D_.begin(), get(graph::vertex_index, graph_)), size_D(size()));
		}


		auto D() const noexcept
		{
			return detail::makeMatrixPtrPropertyMap<OcpVertexDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(D_.begin(), get(graph::vertex_index, graph_)), size_D(size()));
		}


		auto q()
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(q_.begin(), get(graph::vertex_index, graph_)), size_x(size()));
		}


		auto q() const noexcept
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(q_.begin(), get(graph::vertex_index, graph_)), size_x(size()));
		}


		auto r()
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(r_.begin(), get(graph::vertex_index, graph_)), size_u(size()));
		}


		auto r() const noexcept
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(r_.begin(), get(graph::vertex_index, graph_)), size_u(size()));
		}


		auto lx()
		{
			// TODO: check size of put()
			return make_iterator_property_map(lx_.begin(), get(graph::vertex_index, graph_));
		}


		auto lx() const noexcept
		{
			// TODO: check size of put()
			return make_iterator_property_map(lx_.begin(), get(graph::vertex_index, graph_));
		}


		auto ux()
		{
			// TODO: check size of put()
			return make_iterator_property_map(ux_.begin(), get(graph::vertex_index, graph_));
		}


		auto ux() const noexcept
		{
			// TODO: check size of put()
			return make_iterator_property_map(ux_.begin(), get(graph::vertex_index, graph_));
		}


		auto lu()
		{
			// TODO: check size of put()
			return make_iterator_property_map(lu_.begin(), get(graph::vertex_index, graph_));
		}


		auto lu() const noexcept
		{
			// TODO: check size of put()
			return make_iterator_property_map(lu_.begin(), get(graph::vertex_index, graph_));
		}


		auto uu()
		{
			// TODO: check size of put()
			return make_iterator_property_map(uu_.begin(), get(graph::vertex_index, graph_));
		}


		auto uu() const noexcept
		{
			// TODO: check size of put()
			return make_iterator_property_map(uu_.begin(), get(graph::vertex_index, graph_));
		}


		auto ld()
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(lg_.begin(), get(graph::vertex_index, graph_)), size_d(size()));
		}


		auto ld() const noexcept
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(lg_.begin(), get(graph::vertex_index, graph_)), size_d(size()));
		}


		auto ud()
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(ug_.begin(), get(graph::vertex_index, graph_)), size_d(size()));
		}


		auto ud() const noexcept
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(ug_.begin(), get(graph::vertex_index, graph_)), size_d(size()));
		}


		auto A()
		{
			return detail::makeMatrixPtrPropertyMap<OcpEdgeDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(A_.begin(), get(graph::edge_index, graph_)), size_A(size(), graph_));
		}


		auto A() const noexcept
		{
			return detail::makeMatrixPtrPropertyMap<OcpEdgeDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(A_.begin(), get(graph::edge_index, graph_)), size_A(size(), graph_));
		}


		auto B()
		{
			return detail::makeMatrixPtrPropertyMap<OcpEdgeDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(B_.begin(), get(graph::edge_index, graph_)), size_B(size(), graph_));
		}


		auto B() const noexcept
		{
			return detail::makeMatrixPtrPropertyMap<OcpEdgeDescriptor, blaze::CustomMatrix<Real, blaze::unaligned, blaze::unpadded, SO>>(
				make_iterator_property_map(B_.begin(), get(graph::edge_index, graph_)), size_B(size(), graph_));
		}


		auto b()
		{
			return detail::makeVectorPtrPropertyMap<OcpEdgeDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(b_.begin(), get(graph::edge_index, graph_)), size_b(size(), graph_));
		}


		auto b() const noexcept
		{
			return detail::makeVectorPtrPropertyMap<OcpEdgeDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(b_.begin(), get(graph::edge_index, graph_)), size_b(size(), graph_));
		}


		auto x() const noexcept
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(x_.begin(), get(graph::vertex_index, graph_)), size_x(size()));
		}


		auto u() const noexcept
		{
			return detail::makeVectorPtrPropertyMap<OcpVertexDescriptor, blaze::CustomVector<Real, blaze::unaligned, blaze::unpadded>>(
				make_iterator_property_map(u_.begin(), get(graph::vertex_index, graph_)), size_u(size()));
		}

		
	private:
		OcpGraph graph_;
		std::vector<OcpSize> size_;

		std::vector<Real> realPool_;
		std::vector<int> intPool_;

		// --------------------------------
		//
		// HPMPC QP problem data
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
		std::vector<Real *> lb_;
		std::vector<Real *> ub_;
		std::vector<Real *> C_;
		std::vector<Real *> D_;
		std::vector<Real *> lg_;
		std::vector<Real *> ug_;
		std::vector<int> nx_;
		std::vector<int> nu_;
		std::vector<int> nb_;
		std::vector<int> ng_;
		std::vector<int *> hidxb_;

		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lx_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> ux_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> lu_;
		std::vector<blaze::DynamicVector<Real, blaze::columnVector>> uu_;

		// --------------------------------
		//
		// HPMPC QP solution data
		//
		// --------------------------------
		std::vector<Real *> x_;
		std::vector<Real *> u_;
		std::vector<Real *> pi_;
		std::vector<Real *> lam_;
		std::array<Real, 4> infNormRes_;

		/// \brief Number of iterations performed by the QP solver.
		int numIter_ = 0;

		// --------------------------------
		//
		// HPMPC solver data
		//
		// --------------------------------

		// Workspace for HPMPC functions
		std::vector<char> solverWorkspace_;

		// Iteration statistics. HPMPC returns 5 Real numbers per iteration.
		std::vector<IterStat> stat_;

		Real muTol_ = 1e-10;

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

        
            void discover_vertex(OcpVertexDescriptor u, OcpGraph const& g)
            {
				auto const& sz = get(sizeMap_, u);
				auto const max_nb = sz.nu() + sz.nx();

				intCount_ += max_nb;	// idxb
                
				realCount_ += sz.nx() * sz.nx();	// Q
				realCount_ += sz.nu() * sz.nx();	// S
				realCount_ += sz.nu() * sz.nu();	// R 
				realCount_ += sz.nx();	// q
				realCount_ += sz.nu();	// r
				realCount_ += max_nb;	// lb
				realCount_ += max_nb;	// ub

				realCount_ += sz.nc() * sz.nx();	// C
				realCount_ += sz.nc() * sz.nu();	// D
				realCount_ += sz.nc();	// lg
				realCount_ += sz.nc();	// ug

				realCount_ += sz.nx();	// x
				realCount_ += sz.nu();	// u
				realCount_ += 2 * sz.nc() + 2 * max_nb;	// lam
            }


			void tree_edge(OcpEdgeDescriptor e, OcpGraph const& g)
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
            InitVisitor(HpmpcWorkspace& ws)
            :   ws_{ws}
            {
            }

        
            void discover_vertex(OcpVertexDescriptor u, OcpGraph const& g)
            {
				auto const& sz = get(ws_.size(), u);
				auto const max_nb = sz.nu() + sz.nx();

                ws_.nx_.push_back(sz.nx());
				ws_.nu_.push_back(sz.nu());
				ws_.nb_.push_back(max_nb);	// this will be updated during solve()
				ws_.ng_.push_back(sz.nc());
				
				// Allocate idxb array. It will be properly initialized in solve().
				int * const idxb = allocInt(max_nb);
				ws_.hidxb_.push_back(idxb);

				// Init idxb to 0,1,2,... s.t. maximum necessary workspace size is returned.
				for (int i = 0; i < max_nb; ++i)
					idxb[i] = i;

				ws_.Q_.push_back(allocReal(sz.nx() * sz.nx()));
				ws_.S_.push_back(allocReal(sz.nu() * sz.nx()));
				ws_.R_.push_back(allocReal(sz.nu() * sz.nu()));
				ws_.q_.push_back(allocReal(sz.nx()));
				ws_.r_.push_back(allocReal(sz.nu()));

				ws_.lb_.push_back(allocReal(max_nb));	// the array will be updated during solve()
				ws_.lu_.emplace_back(sz.nu(), -inf<Real>());
				ws_.lx_.emplace_back(sz.nx(), -inf<Real>());

				ws_.ub_.push_back(allocReal(max_nb));	// the array will be updated during solve()
				ws_.uu_.emplace_back(sz.nu(), inf<Real>());			
				ws_.ux_.emplace_back(sz.nx(), inf<Real>());

				ws_.C_ .push_back(allocReal(sz.nc() * sz.nx()));
				ws_.D_ .push_back(allocReal(sz.nc() * sz.nu()));
				ws_.lg_.push_back(allocReal(sz.nc()));
				ws_.ug_.push_back(allocReal(sz.nc()));

				ws_.x_.push_back(allocReal(sz.nx()));
				ws_.u_.push_back(allocReal(sz.nu()));
				ws_.lam_.push_back(allocReal(2 * sz.nc() + 2 * max_nb));
            }


			void tree_edge(OcpEdgeDescriptor e, OcpGraph const& g)
			{
				auto const u = source(e, g);
				auto const v = target(e, g);

				if (v != u + 1)
					TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid tree structure in HpmpcWorkspace ctor:	vertices are not sequentially connected"));

				auto const& sz_u = get(ws_.size(), u);
				auto const& sz_v = get(ws_.size(), v);

				ws_.A_.push_back(allocReal(sz_v.nx() * sz_u.nx()));
				ws_.B_.push_back(allocReal(sz_v.nx() * sz_u.nu()));
				ws_.b_.push_back(allocReal(sz_v.nx()));
				ws_.pi_.push_back(allocReal(sz_v.nx()));
			}


			void non_tree_edge(OcpEdgeDescriptor e, OcpGraph const& g)
			{
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid tree structure in HpmpcWorkspace ctor:	non-tree graph structure detected."));
			}

        
        private:
            HpmpcWorkspace& ws_;
			size_t allocatedReal_ = 0;
			size_t allocatedInt_ = 0;


			Real * allocReal(size_t n)
			{
				if (allocatedReal_ + n > ws_.realPool_.size())
					TMPC_THROW_EXCEPTION(std::bad_alloc());

				Real * const ptr = ws_.realPool_.data() + allocatedReal_;
				allocatedReal_ += n;

				return ptr;
			}


			int * allocInt(size_t n)
			{
				if (allocatedInt_ + n > ws_.intPool_.size())
					TMPC_THROW_EXCEPTION(std::bad_alloc());

				int * const ptr = ws_.intPool_.data() + allocatedInt_;
				allocatedInt_ += n;

				return ptr;
			}
        };


		class PopulateBounds
		{
		public:
			PopulateBounds(Real * lb, Real * ub, int * idxb, int& nb)
			:	lb_(lb)
			,	ub_(ub)
			,	idxb_(idxb)
			,	nb_(nb)
			{
				nb_ = 0;
			}


			void operator()(Real l, Real u)
			{
				if (std::isfinite(l) && std::isfinite(u))
				{
					// If both bounds are finite, add i to the bounds index,
					// and copy values to the lb_internal_ and ub_internal_.
					*idxb_++ = count_;
					*lb_++ = l;
					*ub_++ = u;
					++nb_;
				}
				else 
				{
					// Otherwise, check that the values are [-inf, inf]
					if (!(l == -inf<Real>() && u == inf<Real>()))
						TMPC_THROW_EXCEPTION(std::invalid_argument("An invalid QP bound is found. For HPMPC/HPIPM, "
							"the bounds should be either both finite or [-inf, inf]"));
				}

				++count_;
			}


		private:
			Real * lb_;
			Real * ub_;
			int * idxb_;
			int& nb_;
			int count_ = 0;
		};


		// Allocate soverWorkspace_ according to nx, nu, nb, ng etc.
		void allocateSolverWorkspace()
		{
			// Number of QP steps for HPMPC
			auto const N = num_vertices(graph_) - 1;
	
			solverWorkspace_.resize(hpmpc::ip_ocp_hard_tv_work_space_size_bytes<Real>(
				static_cast<int>(N), nx_.data(), nu_.data(), nb_.data(), hidxb_.data(), ng_.data(), static_cast<int>(N)));
		}
	};


	template <typename Real>
    struct RealOf<HpmpcWorkspace<Real>>
    {
        using type = Real;
    };
}
