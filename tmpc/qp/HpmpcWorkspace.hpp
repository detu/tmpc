#pragma once

#include "QpSolverException.hpp"
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpSolutionBase.hpp>
#include <tmpc/qp/OcpQpBase.hpp>
#include <tmpc/qp/QpWorkspaceBase.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

#include "detail/HpxxxVertexData.hpp"
#include "detail/HpxxxEdgeData.hpp"
#include "TreeQpWorkspaceAdaptor.hpp"
//#include "OcpQpVertexElement.hpp"

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/property_map.hpp>

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>
#include <memory>

namespace tmpc
{
	namespace detail
	{
		// Type-parameterized HPMPC interface
		template <typename Real>
		struct Hpmpc
		{
			static int c_order_ip_ocp_hard_tv(
				int *kk, int k_max, Real mu0, Real mu_tol,
				int N, int const *nx, int const *nu, int const *nb, int const * const *hidxb, int const *ng, int N2,
				int warm_start,
				Real const * const *A, Real const * const *B, Real const * const *b,
				Real const * const *Q, Real const * const *S, Real const * const *R, Real const * const *q, Real const * const *r,
				Real const * const *lb, Real const * const *ub,
				Real const * const *C, Real const * const *D, Real const * const *lg, Real const * const *ug,
				Real * const *x, Real * const *u, Real * const *pi, Real * const *lam,
				Real *inf_norm_res,
				void *work0,
				Real *stat);

			static int ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int const * const * hidxb, int const *ng, int N2);
		};
	}


	class HpmpcException 
	:	public std::runtime_error
	{
	public:
		HpmpcException(int code);
		int code() const { return _code;	}

	private:
		int const _code;
	};


	/**
	 * \brief Multistage QP solver using HPMPC
	 *
	 * \tparam <Kernel> the math kernel type
	 */
	template <typename Kernel_>
	class HpmpcWorkspace
	:	public QpWorkspaceBase<HpmpcWorkspace<Kernel_>>
	{
	public:
		using Kernel = Kernel_;
		using SO = StorageOrder::rowMajor;
		using Real = typename Kernel::Real;
		
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
		template <typename InputIterator>
		explicit HpmpcWorkspace(OcpSizeGraph const& graph, unsigned max_iter = 100)
		:	graph_ {graph}
		,	stat_ {max_iter}
		{
			std::fill(infNormRes_.begin(), infNormRes_.end(), sNaN<Real>());

			auto const nv = num_vertices(graph_);
			auto const ne = num_edges(graph_);

			if (nv < 2)
				throw std::invalid_argument("HPMPC needs an at least 2-stages problem");

			if (ne + 1 != num_vertices(graph_))
				throw std::invalid_argument("Number of edges in HPMPC size graph must be 1 less than the number of vertices");

			// Preallocate arrays holding QP vertex data.
			vertexData_.reserve(nv);
	
			nx_.reserve(nv);
			nu_.reserve(nv);
			nb_.reserve(nv);
			ng_.reserve(nv);
			hidxb_.reserve(nv);
	
			A_ .reserve(nv);
			B_ .reserve(nv);
			b_ .reserve(nv);
			Q_ .reserve(nv);
			S_ .reserve(nv);
			R_ .reserve(nv);
			q_ .reserve(nv);
			r_ .reserve(nv);
			lb_.reserve(nv);
			ub_.reserve(nv);
			C_ .reserve(nv);
			D_ .reserve(nv);
			lg_.reserve(nv);
			ug_.reserve(nv);
			
			x_.reserve(nv);
			u_.reserve(nv);
			pi_.reserve(nv);
			lam_.reserve(nv);

			// Preallocate arrays holding QP edge data.
			edgeData_.reserve(ne);
	
			A_ .reserve(ne);
			B_ .reserve(ne);
			b_ .reserve(ne);
			pi_.reserve(ne);

			// Traverse the graph and init vertex and edge data.
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


		void impl_solve()
		{
			// Recalculate bounds indices of each stage, as the bounds might have changed.				
			// Update the nb_ array.
			auto const vd = make_iterator_range(vertices(graph_));
			std::transform(vd.begin(), vd.end(), nb_.begin(), [] (auto& s) -> int
			{
				s.adjustBoundsIndex();
				return s.nb();
			});

			// Number of QP steps for HPMPC
			auto const N = num_vertices(graph_) - 1;

			// What is a good value for mu0?
			Real mu0 = 1.;

			// Call HPMPC
			auto const ret = detail::Hpmpc<Real>::c_order_ip_ocp_hard_tv(&numIter_, maxIter(), mu0, muTol_, N,
					nx_.data(), nu_.data(), nb_.data(), hidxb_.data(), ng_.data(), N, _warmStart ? 1 : 0, A_.data(), B_.data(), b_.data(),
					Q_.data(), S_.data(), R_.data(), q_.data(), r_.data(), lb_.data(), ub_.data(), C_.data(), D_.data(),
					lg_.data(), ug_.data(), x_.data(), u_.data(), pi_.data(), lam_.data(), infNormRes_.data(),
					solverWorkspace_.data(), stat_[0].data());

			if (ret != 0)
				throw HpmpcException(ret);
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
		auto numIter() const 
		{ 
			return numIter_; 
		}


		using VertexPropertyMap = boost::iterator_property_map<
			std::vector<VertexData>::iterator, 


		class VertexPropertyMap
		{
		public:
			VertexPropertyMap(HpmpcWorkspace& ws)
			:	ws_ {ws}
			{				
			}


			VertexData& at(OcpVertexDescriptor k) const
			{
				return ws_.vertexData_.at(k);
			}


		private:
			HpmpcWorkspace& ws_;
		};


		class EdgePropertyMap
		{
		public:
			EdgePropertyMap(HpmpcWorkspace& ws)
			:	ws_ {ws}
			{				
			}


			EdgeData& at(OcpEdgeDescriptor k) const
			{
				return ws_.edgeData_.at(EDGE_INDEX(k));
			}


		private:
			HpmpcWorkspace& ws_;
		};


		template <
			typename EdgeOrVectorPropertyMap,
            Real * (VertexData::* PtrFunc)(),
            std::pair<size_t, size_t> (* SizeFunc)(OcpSizeGraph const&, OcpVertexDescriptor),
			StorageOrder SO = HpmpcWorkspace::SO
        >
		class VertexMatrixPropertyMap
        {
        public:
            VertexMatrixPropertyMap(HpmpcWorkspace& ws)
            :   ws_{ws}
            {
            }


			template <typename T>
            friend void put(VertexMatrixPropertyMap const& pm, OcpVertexDescriptor k, T const& val)
            {
                pm.put(k, val);
            }


            friend decltype(auto) get(VertexMatrixPropertyMap const& pm, OcpVertexDescriptor k)
            {
                return pm.get(k);
            }


        private:
			template <typename T>
            void put(OcpVertexDescriptor k, T const& val) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
				
				CustomMatrix<Real, unaligned, unpadded, SO> lhs(
					(ws_.vertexData_.at(k).*PtrFunc)(), sz.first, sz.second);

                noresize(lhs) = val;
            }


            auto get(OcpVertexDescriptor k) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
                
                CustomMatrix<Real const, unaligned, unpadded, SO> m(
					(ws_.vertexData_.at(k).*PtrFunc)(), sz.first, sz.second);

                return m;
            }


            HpmpcWorkspace& ws_;
        };


		template <
            Real * (VertexData::* PtrFunc)(),
            size_t (* SizeFunc)(OcpSizeGraph const&, OcpVertexDescriptor)
        >
		class VertexVectorPropertyMap
        {
        public:
            VertexVectorPropertyMap(HpmpcWorkspace& ws)
            :   ws_{ws}
            {
            }


			template <typename T>
            friend void put(VertexVectorPropertyMap const& pm, OcpVertexDescriptor k, T const& val)
            {
                pm.put(k, val);
            }


            friend decltype(auto) get(VertexVectorPropertyMap const& pm, OcpVertexDescriptor k)
            {
                return pm.get(k);
            }


        private:
			template <typename T>
            void put(OcpVertexDescriptor k, T const& val) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
				
				CustomVector<Real, unaligned, unpadded, columnVector> lhs(
					(ws_.vertexData_.at(k).*PtrFunc)(), sz);

                noresize(lhs) = val;
            }


            auto get(OcpVertexDescriptor k) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
                
                CustomVector<Real const, unaligned, unpadded, columnVector> val(
					(ws_.vertexData_.at(k).*PtrFunc)(), sz);

                return val;
            }


            HpmpcWorkspace& ws_;
        };


		template <
            Real * (EdgeData::* PtrFunc)(),
            std::pair<size_t, size_t> (* SizeFunc)(OcpSizeGraph const&, OcpEdgeDescriptor),
			StorageOrder SO = HpmpcWorkspace::SO
        >
		class EdgeMatrixPropertyMap
        {
        public:
            EdgeMatrixPropertyMap(HpmpcWorkspace& ws)
            :   ws_{ws}
            {
            }


			template <typename T>
            friend void put(EdgeMatrixPropertyMap const& pm, OcpEdgeDescriptor k, T const& val)
            {
                pm.put(k, val);
            }


            friend decltype(auto) get(EdgeMatrixPropertyMap const& pm, OcpEdgeDescriptor k)
            {
                return pm.get(k);
            }


        private:
			template <typename T>
            void put(OcpEdgeDescriptor k, T const& val) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
				
				CustomMatrix<Real, unaligned, unpadded, SO> lhs(
					(ws_.edgeData_.at(k).*PtrFunc)(), sz.first, sz.second);

                noresize(lhs) = val;
            }


            auto get(OcpEdgeDescriptor k) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
                
                CustomMatrix<Real const, unaligned, unpadded, SO> m(
					(ws_.edgeData_.at(k).*PtrFunc)(), sz.first, sz.second);

                return m;
            }


            HpmpcWorkspace& ws_;
        };


		template <
            Real * (VertexData::* PtrFunc)(),
            size_t (* SizeFunc)(OcpSizeGraph const&, Key)
        >
		class EdgeVectorPropertyMap
        {
        public:
            EdgeVectorPropertyMap(HpmpcWorkspace& ws)
            :   ws_{ws}
            {
            }


			template <typename T>
            friend void put(EdgeVectorPropertyMap const& pm, OcpEdgeDescriptor k, T const& val)
            {
                pm.put(k, val);
            }


            friend decltype(auto) get(EdgeVectorPropertyMap const& pm, OcpEdgeDescriptor k)
            {
                return pm.get(k);
            }


        private:
			template <typename T>
            void put(OcpEdgeDescriptor k, T const& val) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
				
				CustomVector<Real, unaligned, unpadded, columnVector> lhs(
					(ws_.edgeData_.at(k).*PtrFunc)(), sz);

                noresize(lhs) = val;
            }


            auto get(OcpEdgeDescriptor k) const
            {
                auto const sz = SizeFunc(ws_.graph_, k);
                
                CustomVector<Real const, unaligned, unpadded, columnVector> val(
					(ws_.edgeData_.at(k).*PtrFunc)(), sz);

                return val;
            }


            HpmpcWorkspace& ws_;
        };

		
	private:
		using VertexData = detail::HpxxxVertexData<Kernel, SO>;
		using EdgeData = detail::HpxxxEdgeData<Kernel, SO>;

		OcpSizeGraph graph_;
		std::vector<VertexData> vertexData_;
		std::vector<EdgeData> edgeData_;

		// --------------------------------
		//
		// HPMPC QP problem data
		//
		// --------------------------------
		std::vector<Real const *> A_;
		std::vector<Real const *> B_;
		std::vector<Real const *> b_;
		std::vector<Real const *> Q_;
		std::vector<Real const *> S_;
		std::vector<Real const *> R_;
		std::vector<Real const *> q_;
		std::vector<Real const *> r_;
		std::vector<Real const *> lb_;
		std::vector<Real const *> ub_;
		std::vector<Real const *> C_;
		std::vector<Real const *> D_;
		std::vector<Real const *> lg_;
		std::vector<Real const *> ug_;
		std::vector<int> nx_;
		std::vector<int> nu_;
		std::vector<int> nb_;
		std::vector<int> ng_;
		std::vector<int const *> hidxb_;

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


		class InitVisitor
        :   public boost::default_bfs_visitor 
        {
        public:
            InitVisitor(HpmpcWorkspace& ws)
            :   ws_{ws}
            {                
            }

        
            void discover_vertex(OcpVertexDescriptor u, OcpSizeGraph const& g)
            {
				auto const& sz = g[u].size;

                ws_.vertexData_.emplace_back(sz);
				auto& vd = vertexData_.back();

				ws_.nx_.push_back(sz.nx());
				ws_.nu_.push_back(sz.nu());
				ws_.nb_.push_back(vd.nb());
				ws_.ng_.push_back(sz.nc());
				
				ws_.hidxb_.push_back(vd.hidxb_data());

				ws_.Q_.push_back(vd.Q_data());
				ws_.S_.push_back(vd.S_data());
				ws_.R_.push_back(vd.R_data());
				ws_.q_.push_back(vd.q_data());
				ws_.r_.push_back(vd.r_data());
				ws_.lb_.push_back(vd.lb_data());
				ws_.ub_.push_back(vd.ub_data());

				ws_.C_ .push_back(vd.C_data());
				ws_.D_ .push_back(vd.D_data());
				ws_.lg_.push_back(vd.lg_data());
				ws_.ug_.push_back(vd.ug_data());

				ws_.x_.push_back(vd.x_data());
				ws_.u_.push_back(vd.u_data());
				ws_.lam_.push_back(vd.lam_data());
            }


			void tree_edge(OcpEdgeDescriptor e, OcpSizeGraph const& g)
			{
				auto const vertex_from = source(e, g);
				auto const vertex_to = target(e, g);

				if (vertex_to != vertex_from + 1)
					throw std::invalid_argument("Invalid tree structure in HpmpcWorkspace ctor:	vertices are not sequentially connected");

				ws_.edgeData_.emplace_back(g[source(e, g)].size, g[target(e, g)].size);
				auto& ed = vertexData_.back();

				ws_.A_.push_back(ed.A_data());
				ws_.B_.push_back(ed.B_data());
				ws_.b_.push_back(ed.b_data());
				ws_.pi_.push_back(ed.pi_data());
			}


			void non_tree_edge(OcpEdgeDescriptor e, OcpSizeGraph const& g)
			{
				throw std::invalid_argument("Invalid tree structure in HpmpcWorkspace ctor:	non-tree graph structure detected.");
			}

        
        private:
            HpmpcWorkspace& ws_;
        };


		// Allocate soverWorkspace_ according to nx, nu, nb, ng etc.
		void allocateSolverWorkspace()
		{
			// Number of QP steps for HPMPC
			auto const N = num_vertices(graph_) - 1;
	
			solverWorkspace_.resize(detail::Hpmpc<Real>::ip_ocp_hard_tv_work_space_size_bytes(
				static_cast<int>(N), nx_.data(), nu_.data(), nb_.data(), hidxb_.data(), ng_.data(), static_cast<int>(N)));
		}
	};
}
