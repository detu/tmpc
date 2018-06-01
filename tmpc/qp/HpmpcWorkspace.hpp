#pragma once

#include "QpSolverException.hpp"
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpSolutionBase.hpp>
#include <tmpc/qp/OcpQpBase.hpp>
#include <tmpc/qp/QpWorkspaceBase.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

#include "detail/HpxxxStage.hpp"

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>
#include <memory>

namespace tmpc
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
		using Real = typename Kernel::Real;
		
		// Iteration statistics. HPMPC returns 5 Real numbers per iteration:
		// step length for predictor and corrector, 
		// centering parameter, duality measure for predictor and corrector.
		using IterStat = std::array<Real, 5>;

		using Stage = detail::HpxxxStage<Kernel, StorageOrder::rowMajor>;

		class ProblemIterator
		:	public boost::iterator_adaptor<
				ProblemIterator	// derived
			,	typename std::vector<Stage>::iterator	// base
			,	OcpQpBase<Stage>&	// value
			,	boost::random_access_traversal_tag	// category of traversal
			,	OcpQpBase<Stage>&	// reference
			>
		{
		public:
			ProblemIterator() = default;
			
			ProblemIterator(typename ProblemIterator::base_type const& p)
			:	ProblemIterator::iterator_adaptor_(p)
			{			
			}
		};

		class ConstProblemIterator
		:	public boost::iterator_adaptor<
				ConstProblemIterator	// derived
			,	typename std::vector<Stage>::const_iterator	// base
			,	OcpQpBase<Stage> const&	// value
			,	boost::random_access_traversal_tag	// category of traversal
			,	OcpQpBase<Stage> const&	// reference
			>
		{
		public:
			ConstProblemIterator() = default;
			
			ConstProblemIterator(typename ConstProblemIterator::base_type const& p)
			:	ConstProblemIterator::iterator_adaptor_(p)
			{			
			}
		};
	
		class ConstSolutionIterator
		:	public boost::iterator_adaptor<
				ConstSolutionIterator	// derived
			,	typename std::vector<Stage>::const_iterator	// base
			,	OcpSolutionBase<Stage> const&	// value
			,	boost::random_access_traversal_tag	// category of traversal
			,	OcpSolutionBase<Stage> const&	// reference
			>
		{
		public:
			ConstSolutionIterator() = default;
			
			ConstSolutionIterator(typename ConstSolutionIterator::base_type const& p)
			:	ConstSolutionIterator::iterator_adaptor_(p)
			{			
			}
		};


		std::string impl_solverName() const
		{
			return "HPMPC";
		}

	
		boost::iterator_range<ProblemIterator> problem()
		{
			return boost::iterator_range<ProblemIterator>(stage_.begin(), stage_.end());
		}
	
		boost::iterator_range<ConstProblemIterator> problem() const
		{
			return boost::iterator_range<ConstProblemIterator>(stage_.begin(), stage_.end());
		}
	
		
		boost::iterator_range<ConstSolutionIterator> impl_solution() const
		{
			return boost::iterator_range<ConstSolutionIterator>(stage_.begin(), stage_.end());
		}

		boost::iterator_range<typename std::vector<IterStat>::const_iterator> stat() const
		{
			return boost::make_iterator_range(stat_.begin(), stat_.begin() + numIter_);
		}

		/**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <typename InputIterator>
		explicit HpmpcWorkspace(InputIterator size_first, InputIterator size_last, unsigned max_iter = 100)
		:	stat_ {max_iter}
		{
			std::fill(infNormRes_.begin(), infNormRes_.end(), sNaN<Real>());

			auto const n_stages = std::distance(size_first, size_last);

			if (n_stages < 2)
				throw std::invalid_argument("HPMPC needs at least 2 stages problem");

			preallocateStages(n_stages);

			for (auto sz = size_first; sz != size_last; ++sz)
				addStage(*sz, sz + 1 != size_last ? sz[1].nx() : 0);

			allocateSolverWorkspace();
		}

		template <typename IteratorRange>
		explicit HpmpcWorkspace(IteratorRange sz, int max_iter = 100)
		:	HpmpcWorkspace(sz.begin(), sz.end(), max_iter)
		{
		}

		explicit HpmpcWorkspace(std::initializer_list<OcpSize> sz, int max_iter = 100)
		:	HpmpcWorkspace(sz.begin(), sz.end(), max_iter)
		{
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
			if (stage_.size() > 0)
			{
				// Recalculate bounds indices of each stage, as the bounds might have changed.				
				// Update the nb_ array.
				std::transform(stage_.begin(), stage_.end(), nb_.begin(), [] (Stage& s) -> int
				{
					s.adjustBoundsIndex();
					return s.nb();
				});
	
				// Number of QP steps for HPMPC
				auto const N = stage_.size() - 1;
	
				// What is a good value for mu0?
				Real mu0 = 1.;
	
				// Call HPMPC
				auto const ret = Hpmpc<Real>::c_order_ip_ocp_hard_tv(&numIter_, maxIter(), mu0, muTol_, N,
						nx_.data(), nu_.data(), nb_.data(), hidxb_.data(), ng_.data(), N, _warmStart ? 1 : 0, A_.data(), B_.data(), b_.data(),
						Q_.data(), S_.data(), R_.data(), q_.data(), r_.data(), lb_.data(), ub_.data(), C_.data(), D_.data(),
						lg_.data(), ug_.data(), x_.data(), u_.data(), pi_.data(), lam_.data(), infNormRes_.data(),
						solverWorkspace_.data(), stat_[0].data());
	
				if (ret != 0)
				{
					throw HpmpcException(ret);
				}
			}
		}

		std::size_t maxIter() const noexcept { return stat_.size(); }

		/// \brief Set max number of iterations
		void maxIter(size_t val)
		{
			stat_.resize(val);
		}

		Real muTol() const noexcept { return muTol_; }

		void muTol(Real val)
		{
			if (val <= 0.)
				throw std::invalid_argument("mu tolerance for hpmpc must be positive");

			muTol_ = val;
		}

		/// \brief Get number of iterations performed by the QP solver.
		unsigned numIter() const { return numIter_; }
		
	private:
		// Stores stage data
		std::vector<Stage> stage_;

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

		// Preallocate arrays holding QP stage data.
		void preallocateStages(size_t nt)
		{
			stage_.reserve(nt);
	
			nx_.reserve(nt);
			nu_.reserve(nt);
			nb_.reserve(nt);
			ng_.reserve(nt);
			hidxb_.reserve(nt);
	
			A_ .reserve(nt);
			B_ .reserve(nt);
			b_ .reserve(nt);
			Q_ .reserve(nt);
			S_ .reserve(nt);
			R_ .reserve(nt);
			q_ .reserve(nt);
			r_ .reserve(nt);
			lb_.reserve(nt);
			ub_.reserve(nt);
			C_ .reserve(nt);
			D_ .reserve(nt);
			lg_.reserve(nt);
			ug_.reserve(nt);
			
			x_.reserve(nt);
			u_.reserve(nt);
			pi_.reserve(nt);
			lam_.reserve(nt);
		}

		// Add one stage with specified sizes at the end of the stage sequence.
		void addStage(OcpSize const& sz, size_t nx_next)
		{
			stage_.emplace_back(sz, nx_next);
			auto& st = stage_.back();

			nx_.push_back(sz.nx());
			nu_.push_back(sz.nu());
			nb_.push_back(st.nb());
			ng_.push_back(sz.nc());
			
			hidxb_.push_back(st.hidxb_data());

			A_.push_back(st.A_data());
			B_.push_back(st.B_data());
			b_.push_back(st.b_data());

			Q_.push_back(st.Q_data());
			S_.push_back(st.S_data());
			R_.push_back(st.R_data());
			q_.push_back(st.q_data());
			r_.push_back(st.r_data());
			lb_.push_back(st.lb_data());
			ub_.push_back(st.ub_data());

			C_ .push_back(st.C_data());
			D_ .push_back(st.D_data());
			lg_.push_back(st.lg_data());
			ug_.push_back(st.ug_data());

			x_.push_back(st.x_data());
			u_.push_back(st.u_data());
			pi_.push_back(st.pi_data());
			lam_.push_back(st.lam_data());
		}

		// Allocate soverWorkspace_ according to nx, nu, nb, ng etc.
		void allocateSolverWorkspace()
		{
			// Number of QP steps for HPMPC
			auto const N = stage_.size() - 1;
	
			solverWorkspace_.resize(Hpmpc<Real>::ip_ocp_hard_tv_work_space_size_bytes(
				static_cast<int>(N), nx_.data(), nu_.data(), nb_.data(), hidxb_.data(), ng_.data(), static_cast<int>(N)));
		}
	};
}
