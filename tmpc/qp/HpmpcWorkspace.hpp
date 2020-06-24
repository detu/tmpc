#pragma once

#include <tmpc/qp/QpSolverException.hpp>
#include <tmpc/qp/HpmpcErrorInfo.hpp>
#include <tmpc/qp/HpmpcIterationInfo.hpp>
#include <tmpc/qp/HpmpcResidualNorm.hpp>

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSolutionBase.hpp>
#include <tmpc/qp/OcpQpBase.hpp>
//#include <tmpc/qp/QpWorkspaceBase.hpp>
#include <tmpc/qp/hpmpc/Hpmpc.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/math/UnpaddedMatrix.hpp>

//#include "detail/HpxxxVertexData.hpp"
//#include "detail/HpxxxEdgeData.hpp"
#include <tmpc/property_map/BundlePropertyMap.hpp>
//#include "TreeQpWorkspaceAdaptor.hpp"
//#include "OcpQpVertexElement.hpp"
#include <tmpc/property_map/PropertyMap.hpp>
#include <tmpc/Exception.hpp>
#include <tmpc/Traits.hpp>

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <vector>
#include <array>
#include <memory>


namespace tmpc
{
	/**
	 * \brief Multistage QP solver using HPMPC
	 *
	 * \tparam <Kernel> the math kernel type
	 */
	template <typename Real>
	class HpmpcWorkspace
	{
	public:
		static auto constexpr SO = blaze::rowMajor;
		
		
		std::string impl_solverName() const
		{
			return "HPMPC";
		}

	
		auto const& iterationInfo() const noexcept
		{
			return iterationInfo_;
		}


        auto const& residualNorm() const noexcept
        {
            return infNormRes_;
        }


		/**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <typename SizeMap>
		explicit HpmpcWorkspace(OcpGraph const& graph, SizeMap size, size_t max_iter = 100)
		:	graph_ {graph}
		,	size_(num_vertices(graph))
        ,	maxIter_(max_iter)
		{
			copyProperty(size, iterator_property_map(size_.begin(), get(graph::vertex_index, graph_)), graph::vertices(graph_));

            // Init vertex and edge data
			initData();

            // Initialize hpmpc data structures
            initHpmpcData();
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
			// Number of QP steps for HPMPC
			auto const N = num_vertices(graph_) - 1;

			// What is a good value for mu0?
			Real mu0 = 1.;

			// Call HPMPC
			try
			{
				hpmpc::c_order_ip_ocp_hard_tv(numIter_, maxIter_, mu0, muTol_, N,
					nx_.data(), nu_.data(), nb_.data(), hidxb_.data(), ng_.data(), N, _warmStart ? 1 : 0, A_.data(), B_.data(), b_.data(),
					Q_.data(), S_.data(), R_.data(), q_.data(), r_.data(), lb_.data(), ub_.data(), C_.data(), D_.data(),
					lg_.data(), ug_.data(), x_.data(), u_.data(), pi_.data(), lam_.data(), infNormRes_.data(),
					solverWorkspace_.get(), reinterpret_cast<Real *>(iterationInfo_.data()));
			}
			catch (QpSolverException const&)
			{
				// Resize iterationInfo_ to correct size even if the QP solver has failed,
				// because the iteration info has been filled.
				iterationInfo_.resize(numIter_);
				throw;
			}

			iterationInfo_.resize(numIter_);
		}


		std::size_t maxIter() const noexcept 
		{ 
			return maxIter_; 
		}


		/// \brief Set max number of iterations
		void maxIter(size_t val)
		{
			iterationInfo_.reserve(val);
			maxIter_ = val;
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
            return BundlePropertyMap(&VertexData::Q_, vertexProperties());
        }


        auto Q() const
        {
            return BundlePropertyMap(&VertexData::Q_, vertexProperties());
        }


        auto R()
        {
            return BundlePropertyMap(&VertexData::R_, vertexProperties());
        }


        auto R() const
        {
            return BundlePropertyMap(&VertexData::R_, vertexProperties());
        }


        auto S()
        {
            return BundlePropertyMap(&VertexData::S_, vertexProperties());
        }


        auto S() const
        {
            return BundlePropertyMap(&VertexData::S_, vertexProperties());
        }


        auto q()
        {
            return BundlePropertyMap(&VertexData::q_, vertexProperties());
        }


        auto q() const
        {
            return BundlePropertyMap(&VertexData::q_, vertexProperties());
        }


        auto r()
        {
            return BundlePropertyMap(&VertexData::r_, vertexProperties());
        }


        auto r() const
        {
            return BundlePropertyMap(&VertexData::r_, vertexProperties());
        }


        auto lx()
        {
            return BundlePropertyMap(&VertexData::lx_, vertexProperties());
        }


        auto lx() const
        {
            return BundlePropertyMap(&VertexData::lx_, vertexProperties());
        }


        auto ux()
        {
            return BundlePropertyMap(&VertexData::ux_, vertexProperties());
        }


        auto ux() const
        {
            return BundlePropertyMap(&VertexData::ux_, vertexProperties());
        }


        auto lu()
        {
            return BundlePropertyMap(&VertexData::lu_, vertexProperties());
        }


        auto lu() const
        {
            return BundlePropertyMap(&VertexData::lu_, vertexProperties());
        }


        auto uu()
        {
            return BundlePropertyMap(&VertexData::uu_, vertexProperties());
        }


        auto uu() const
        {
            return BundlePropertyMap(&VertexData::uu_, vertexProperties());
        }


        auto C()
        {
            return BundlePropertyMap(&VertexData::C_, vertexProperties());
        }


        auto C() const
        {
            return BundlePropertyMap(&VertexData::C_, vertexProperties());
        }


        auto D()
        {
            return BundlePropertyMap(&VertexData::D_, vertexProperties());
        }


        auto D() const
        {
            return BundlePropertyMap(&VertexData::D_, vertexProperties());
        }


        auto ld()
        {
            return BundlePropertyMap(&VertexData::lg_, vertexProperties());
        }


        auto ld() const
        {
            return BundlePropertyMap(&VertexData::lg_, vertexProperties());
        }


        auto ud()
        {
            return BundlePropertyMap(&VertexData::ug_, vertexProperties());
        }


        auto ud() const
        {
            return BundlePropertyMap(&VertexData::ug_, vertexProperties());
        }


		auto A()
        {
            return BundlePropertyMap(&EdgeData::A_, edgeProperties());
        }


        auto A() const
        {
            return BundlePropertyMap(&EdgeData::A_, edgeProperties());
        }


        auto B()
        {
            return BundlePropertyMap(&EdgeData::B_, edgeProperties());
        }


        auto B() const
        {
            return BundlePropertyMap(&EdgeData::B_, edgeProperties());
        }


        auto b()
        {
            return BundlePropertyMap(&EdgeData::b_, edgeProperties());
        }


        auto b() const
        {
            return BundlePropertyMap(&EdgeData::b_, edgeProperties());
        }


        auto x()
        {
            return BundlePropertyMap(&VertexData::x_, vertexProperties());
        }


        auto x() const
        {
            return BundlePropertyMap(&VertexData::x_, vertexProperties());
        }


        auto u()
        {
            return BundlePropertyMap(&VertexData::u_, vertexProperties());
        }


        auto u() const
        {
            return BundlePropertyMap(&VertexData::u_, vertexProperties());
        }


        auto lam_lx()
        {
            return BundlePropertyMap(&VertexData::lam_lx_, vertexProperties());
        }


        auto lam_lx() const
        {
            return BundlePropertyMap(&VertexData::lam_lx_, vertexProperties());
        }


        auto lam_ux()
        {
            return BundlePropertyMap(&VertexData::lam_ux_, vertexProperties());
        }


        auto lam_ux() const
        {
            return BundlePropertyMap(&VertexData::lam_ux_, vertexProperties());
        }


        auto lam_lu()
        {
            return BundlePropertyMap(&VertexData::lam_lu_, vertexProperties());
        }


        auto lam_lu() const
        {
            return BundlePropertyMap(&VertexData::lam_lu_, vertexProperties());
        }


        auto lam_uu()
        {
            return BundlePropertyMap(&VertexData::lam_uu_, vertexProperties());
        }


        auto lam_uu() const
        {
            return BundlePropertyMap(&VertexData::lam_uu_, vertexProperties());
        }


        auto lam_ld()
        {
            return BundlePropertyMap(&VertexData::lam_ld_, vertexProperties());
        }


        auto lam_ld() const
        {
            return BundlePropertyMap(&VertexData::lam_ld_, vertexProperties());
        }


        auto lam_ud()
        {
            return BundlePropertyMap(&VertexData::lam_ud_, vertexProperties());
        }


        auto lam_ud() const
        {
            return BundlePropertyMap(&VertexData::lam_ud_, vertexProperties());
        }


        auto pi()
        {
            return BundlePropertyMap(&EdgeData::pi_, edgeProperties());
        }


        auto pi() const
        {
            return BundlePropertyMap(&EdgeData::pi_, edgeProperties());
        }

		
	private:
		struct VertexData
        {
            VertexData(OcpSize const& sz)
            :   size_ {sz}
            ,   Q_(sz.nx(), sz.nx())
            ,   R_(sz.nu(), sz.nu())
            ,   S_(sz.nu(), sz.nx())
            ,   q_(sz.nx())
            ,   r_(sz.nu())
			,	lb_(sz.nx() + sz.nu())
			,	ub_(sz.nx() + sz.nu())
            ,   lx_(subvector(lb_, sz.nu(), sz.nx()))
            ,   ux_(subvector(ub_, sz.nu(), sz.nx()))
            ,   lu_(subvector(lb_, 0, sz.nu()))
            ,   uu_(subvector(ub_, 0, sz.nu()))
            ,   C_(sz.nc(), sz.nx())
            ,   D_(sz.nc(), sz.nu())
            ,   lg_(sz.nc())
            ,   ug_(sz.nc())
            ,   x_(blaze::ZeroVector<Real>(sz.nx()))
            ,   u_(blaze::ZeroVector<Real>(sz.nu()))
			,	lam_(sz.nu() + sz.nx() + sz.nc() + sz.nu() + sz.nx() + sz.nc())
			//	I am not sure about the order of lambdas returned by HPMPC
            ,   lam_lu_(subvector(lam_, 0, sz.nu()))
            ,   lam_lx_(subvector(lam_, sz.nu(), sz.nx()))
            ,   lam_ld_(subvector(lam_, sz.nu() + sz.nx(), sz.nc()))
            ,   lam_uu_(subvector(lam_, sz.nu() + sz.nx() + sz.nc(), sz.nu()))
            ,   lam_ux_(subvector(lam_, sz.nu() + sz.nx() + sz.nc() + sz.nu(), sz.nx()))
            ,   lam_ud_(subvector(lam_, sz.nu() + sz.nx() + sz.nc() + sz.nu() + sz.nx(), sz.nc()))
			,	idxb_ {new int[sz.nx() + sz.nu()]}
            {
				// Init idxb to 0,1,2,...
				for (int i = 0; i < sz.nx() + sz.nu(); ++i)
					idxb_[i] = i;
            }


            // We don't want VertexData to be copied.
            VertexData(VertexData const&) = delete;
            

            // We don't want VertexData to be moved, but std::vector::reserve() requires move ctor.
            // Although never used, it needs to be defined.
            VertexData(VertexData&& rhs)
            :   size_ {rhs.size_}
            ,   Q_(std::move(rhs.Q_))
            ,   R_(std::move(rhs.R_))
            ,   S_(std::move(rhs.S_))
            ,   q_(std::move(rhs.q_))
            ,   r_(std::move(rhs.r_))
			,	lb_(std::move(rhs.lb_))
			,	ub_(std::move(rhs.ub_))
            ,   lx_(subvector(lb_, size_.nu(), size_.nx()))
            ,   ux_(subvector(ub_, size_.nu(), size_.nx()))
            ,   lu_(subvector(lb_, 0, size_.nu()))
            ,   uu_(subvector(ub_, 0, size_.nu()))
            ,   C_(std::move(rhs.C_))
            ,   D_(std::move(rhs.D_))
            ,   lg_(std::move(rhs.lg_))
            ,   ug_(std::move(rhs.ug_))
            ,   x_(std::move(rhs.x_))
            ,   u_(std::move(rhs.u_))
			,	lam_(std::move(rhs.lam_))
			//	I am not sure about the order of lambdas returned by HPMPC
            ,   lam_lu_(subvector(lam_, 0, size_.nu()))
            ,   lam_lx_(subvector(lam_, size_.nu(), size_.nx()))
            ,   lam_ld_(subvector(lam_, size_.nu() + size_.nx(), size_.nc()))
            ,   lam_uu_(subvector(lam_, size_.nu() + size_.nx() + size_.nc(), size_.nu()))
            ,   lam_ux_(subvector(lam_, size_.nu() + size_.nx() + size_.nc() + size_.nu(), size_.nx()))
            ,   lam_ud_(subvector(lam_, size_.nu() + size_.nx() + size_.nc() + size_.nu() + size_.nx(), size_.nc()))
			,	idxb_ {std::move(rhs.idxb_)}
            {
                TMPC_THROW_EXCEPTION(std::logic_error(
                    "If we get here, probably vertexData_ reallocation has occured -- this should never happen!"));
            }


            OcpSize const size_;

            UnpaddedMatrix<Real, SO> Q_;
            UnpaddedMatrix<Real, SO> R_;
            UnpaddedMatrix<Real, SO> S_;
            blaze::DynamicVector<Real> q_;
            blaze::DynamicVector<Real> r_;

			blaze::DynamicVector<Real> lb_;
            blaze::DynamicVector<Real> ub_;

			using Subvector = decltype(blaze::subvector(lb_, 0, 1));

            Subvector lx_;
            Subvector ux_;
            Subvector lu_;
            Subvector uu_;

            UnpaddedMatrix<Real, SO> C_;
            UnpaddedMatrix<Real, SO> D_;
            blaze::DynamicVector<Real> lg_;
            blaze::DynamicVector<Real> ug_;

            blaze::DynamicVector<Real> x_;
            blaze::DynamicVector<Real> u_;
			blaze::DynamicVector<Real> lam_;
            Subvector lam_lu_;
            Subvector lam_lx_;
            Subvector lam_ld_;
            Subvector lam_uu_;
            Subvector lam_ux_;
            Subvector lam_ud_;

			std::unique_ptr<int []> idxb_;
        };


        struct EdgeData
        {
            EdgeData(OcpSize const& size_src, OcpSize const& size_dst)
            :   A_(size_dst.nx(), size_src.nx())
            ,   B_(size_dst.nx(), size_src.nu())
			,   b_(size_dst.nx())
            ,   pi_(size_dst.nx())
            {
            }

            // We don't want EdgeData to be copied.
            EdgeData(EdgeData const&) = delete;

            
            // We don't want EdgeData to be moved, but std::vector::reserve() requires move ctor.
            // Although never used, it needs to be defined.
            EdgeData(EdgeData&& rhs)
            :   A_(std::move(rhs.A_))
            ,   B_(std::move(rhs.B_))
			,   b_(std::move(rhs.b_))
            ,   pi_(std::move(rhs.pi_))
            {
                TMPC_THROW_EXCEPTION(std::logic_error(
                    "If we get here, probably edgeData_ reallocation has occured -- this should never happen!"));
            };


            UnpaddedMatrix<Real, SO> A_;
            UnpaddedMatrix<Real, SO> B_;
			
			/*
            UnpaddedMatrix<Real, SO> BA_;

            using Submatrix = decltype(blaze::submatrix(BA_, 0, 0, 1, 1));

            Submatrix A_;
            Submatrix B_;
			*/
			
            blaze::DynamicVector<Real> b_;
            blaze::DynamicVector<Real> pi_;
        };


		auto vertexProperties()
        {
            return make_iterator_property_map(begin(vertexData_), get(graph::vertex_index, graph_));
        }


        auto vertexProperties() const
        {
            return make_iterator_property_map(begin(vertexData_), get(graph::vertex_index, graph_));
        }


        auto edgeProperties()
        {
            return make_iterator_property_map(begin(edgeData_), get(graph::edge_index, graph_));
        }


        auto edgeProperties() const
        {
            return make_iterator_property_map(begin(edgeData_), get(graph::edge_index, graph_));
        }


        void initData()
        {
            // Preallocate arrays holding QP edge and vertex data.
            vertexData_.reserve(num_vertices(graph_));
            edgeData_.reserve(num_edges(graph_));
			
			// Traverse the graph and init vertex and edge data pointers.
			breadth_first_search(graph_, vertex(0, graph_), visitor(InitVisitor(*this)));
        }


        /// @brief Preallocate arrays holding HPMPC data pointers.
        void initHpmpcData()
        {
            // Init vertex data
            auto const nv = vertexData_.size();
            if (nv < 2)
				TMPC_THROW_EXCEPTION(std::invalid_argument("HPMPC needs a graph with at least 2 nodes"));

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

			for (auto v : graph::vertices(graph_))
            {
                auto const& sz = get(size(), v);
                auto& vd = get(vertexProperties(), v);

                nx_.push_back(sz.nx());
				nu_.push_back(sz.nu());
				nb_.push_back(sz.nu() + sz.nx());
				ng_.push_back(sz.nc());

				hidxb_.push_back(vd.idxb_.get());

				Q_.push_back(vd.Q_.data());
				S_.push_back(vd.S_.data());
				R_.push_back(vd.R_.data());
				q_.push_back(vd.q_.data());
				r_.push_back(vd.r_.data());

				lb_.push_back(vd.lb_.data());
				ub_.push_back(vd.ub_.data());

				C_ .push_back(vd.C_.data());
				D_ .push_back(vd.D_.data());
				lg_.push_back(vd.lg_.data());
				ug_.push_back(vd.ug_.data());

				x_.push_back(vd.x_.data());
				u_.push_back(vd.u_.data());
				lam_.push_back(vd.lam_.data());
            }


            // Init edge data
            auto const ne = edgeData_.size();
            if (ne + 1 != num_vertices(graph_))
				TMPC_THROW_EXCEPTION(std::invalid_argument("Number of edges in HPMPC size graph must be 1 less than the number of vertices"));

            A_.reserve(ne);
			B_.reserve(ne);
			b_.reserve(ne);
			pi_.reserve(ne);

            for (auto e : graph::edges(graph_))
            {
                auto& ed = get(edgeProperties(), e);

				A_.push_back(ed.A_.data());
				B_.push_back(ed.B_.data());
				b_.push_back(ed.b_.data());
				pi_.push_back(ed.pi_.data());
            }


            // Allocate HPMPC working memory according to nx, nu, nb, ng etc.
			// Number of QP steps for HPMPC
			auto const N = nv - 1;
			solverWorkspace_.reset(new char[hpmpc::ip_ocp_hard_tv_work_space_size_bytes<Real>(
				N, nx_.data(), nu_.data(), nb_.data(), hidxb_.data(), ng_.data(), N)]);

            iterationInfo_.reserve(maxIter_);            
			std::fill(infNormRes_.begin(), infNormRes_.end(), sNaN<Real>());
        }


		OcpGraph graph_;
		std::vector<OcpSize> size_;

		std::vector<VertexData> vertexData_;
		std::vector<EdgeData> edgeData_;

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
		HpmpcResidualNorm<Real> infNormRes_;

		/// \brief Number of iterations performed by the QP solver.
		int numIter_ = 0;

		// --------------------------------
		//
		// HPMPC solver data
		//
		// --------------------------------

		// Workspace for HPMPC functions
		std::unique_ptr<char []> solverWorkspace_;

		// Maximum allowed number of iterations
		size_t maxIter_;

		// Iteration info.
		std::vector<HpmpcIterationInfo<Real>> iterationInfo_;

		Real muTol_ = 1e-10;

		// Warmstarting disabled on purpose.
		// On AMD K8 (hpmpc compiled for SSE3), WITHOUT warmstarting it is significantly
		// FASTER (9ms vs 14ms per time step) than with warmstarting. I am curious why.
		bool _warmStart = false;


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
				ws_.vertexData_.emplace_back(get(ws_.size(), u));
            }


			void tree_edge(OcpEdgeDescriptor e, OcpGraph const& g)
			{
				auto const u = source(e, g);
				auto const v = target(e, g);

				if (v != u + 1)
					TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid tree structure in HpmpcWorkspace ctor: vertices are not sequentially connected."));

				auto const& sz_u = get(ws_.size(), u);
				auto const& sz_v = get(ws_.size(), v);
				ws_.edgeData_.emplace_back(sz_u, sz_v);
			}


			void non_tree_edge(OcpEdgeDescriptor e, OcpGraph const& g)
			{
				TMPC_THROW_EXCEPTION(std::invalid_argument("Invalid tree structure in HpmpcWorkspace ctor: non-tree graph structure detected."));
			}

        
        private:
            HpmpcWorkspace& ws_;
        };
	};


	template <typename Real>
    struct RealOf<HpmpcWorkspace<Real>>
    {
        using type = Real;
    };
}
