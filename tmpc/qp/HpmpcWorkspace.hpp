#pragma once

#include "QpSolverException.hpp"
#include "QpSize.hpp"
#include "QpStageSolutionBase.hpp"
#include "QpStageBase.hpp"

#include <tmpc/Math.hpp>
#include <tmpc/matrix/StorageOrder.hpp>
#include <tmpc/matrix/TransposeFlag.hpp>

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <array>
#include <cmath>

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
	:	public QpSolverException
	{
	public:
		HpmpcException(int code);

		int code() const { return _code;	}
		char const * what() const noexcept override { return msg_.c_str(); }

	private:
		int const _code;
		std::string const msg_;
	};

	/**
	 * \brief Multistage QP solver using HPMPC
	 *
	 * \tparam <Kernel> the math kernel type
	 */
	template <typename Kernel_>
	class HpmpcWorkspace
	{
	public:
		using Kernel = Kernel_;
		using Real = typename Kernel::Real;
		
		// Iteration statistics. HPMPC returns 5 Real numbers per iteration:
		// step length for predictor and corrector, 
		// centering parameter, duality measure for predictor and corrector.
		using IterStat = std::array<Real, 5>;

		class Stage
		:	public QpStageSolutionBase<Stage>
		,	public QpStageBase<Stage>
		{
		public:
			static auto constexpr storageOrder = rowMajor;
			using Matrix = typename Kernel::DynamicMatrix;
			using Vector = typename Kernel::DynamicVector;

			Stage(QpSize const& sz, size_t nx_next)
			:	size_(sz)
			,	hidxb_(sz.nu() + sz.nx())
			// Initialize all numeric data to NaN so that if an uninitialized object
			// by mistake used in calculations is easier to detect.
			,	Q_(sz.nx(), sz.nx(), sNaN<Real>())
			,	R_(sz.nu(), sz.nu(), sNaN<Real>())
			,	S_(sz.nu(), sz.nx(), sNaN<Real>())	// <-- HPMPC convention for S is [nu, nx] (the corresponding cost term is u' * S_{hpmpc} * x)
			,	q_(sz.nx(), sNaN<Real>())
			,	r_(sz.nu(), sNaN<Real>())
			,	A_(nx_next, sz.nx(), sNaN<Real>())
			,	B_(nx_next, sz.nu(), sNaN<Real>())
			,	b_(nx_next, sNaN<Real>())
			,	C_(sz.nc(), sz.nx(), sNaN<Real>())
			,	D_(sz.nc(), sz.nu(), sNaN<Real>())
			,	lb_(sz.nu() + sz.nx(), sNaN<Real>())
			,	ub_(sz.nu() + sz.nx(), sNaN<Real>())
			,	lb_internal_(sz.nu() + sz.nx(), sNaN<Real>())
			,	ub_internal_(sz.nu() + sz.nx(), sNaN<Real>())
			,	lbd_(sz.nc(), sNaN<Real>())
			,	ubd_(sz.nc(), sNaN<Real>())
			,	x_(sz.nx(), sNaN<Real>())
			,	u_(sz.nu(), sNaN<Real>())
			,	pi_(nx_next, sNaN<Real>())
			,	lam_(2 * sz.nc() + 2 * (sz.nx() + sz.nu()), sNaN<Real>())
			{
				// hidxb is initialized to its maximum size, s.t. nb == nx + nu.
				// This is necessary so that the solver workspace memory is calculated as its maximum when allocated.
				int n = 0;
				std::generate(hidxb_.begin(), hidxb_.end(), [&n] { return n++; });
			}

			Stage(Stage const&) = default;
			Stage(Stage &&) = default;

			template <typename Expr>
			Stage& operator=(Expr const& rhs)
			{
				assign(*this, rhs);
				return *this;
			}

			auto const& A() const { return A_; }
			template <typename T> void A(const T& a) { full(A_) = a; }

			auto const& B() const { return B_; }
			template <typename T> void B(const T& b) { full(B_) = b; }

			auto const& b() const { return b_; }
			template <typename T> void b(const T& b) { full(b_) = b; }

			auto const& C() const { return C_; }
			template <typename T> void C(const T& c) { full(C_) = c; }

			auto const& D() const { return D_; }
			template <typename T> void D(const T& d) { full(D_) = d; }

			auto const& lbd() const {	return lbd_; }
			template <typename T> void lbd(const T& lbd) { full(lbd_) = lbd; }

			auto lbu() const { return subvector(lb_, 0, size_.nu());	}		
			
			// TODO: consider setting both upper and lower bounds at the same time.
			// Maybe create a Bounds class?
			template <typename T> void lbu(const T& lbu) 
			{ 
				subvector(lb_, 0, size_.nu()) = lbu;
			}

			auto lbx() const { return subvector(lb_, size_.nu(), size_.nx()); }
			template <typename T> void lbx(const T& lbx) 
			{ 
				subvector(lb_, size_.nu(), size_.nx()) = lbx; 
			}

			auto const& Q() const { return Q_; }
			template <typename T> void Q(const T& q) { full(Q_) = q; }

			auto const& R() const { return R_; }
			template <typename T> void R(const T& r) { full(R_) = r; }

			// HPMPC convention for S is [nu, nx], therefore the trans().
			auto S() const { return trans(S_); }
			void S(Real v) { S_ = v; }
			template <typename T> void S(const T& s) { full(S_) = trans(s); }

			auto const& q() const { return q_; }
			template <typename T> void q(const T& q) { full(q_) = q; }

			auto const& r() const { return r_; }
			template <typename T> void r(const T& r) { full(r_) = r; }

			auto const& ubd() const { return ubd_; }
			template <typename T> void ubd(const T& ubd) { full(ubd_) = ubd; }

			auto ubu() const { return subvector(ub_, 0, size_.nu()); }
			template <typename T> void ubu(const T& ubu) 
			{ 
				subvector(ub_, 0, size_.nu()) = ubu;
			}

			auto ubx() const { return subvector(ub_, size_.nu(), size_.nx()); }
			template <typename T> void ubx(const T& ubx) 
			{ 
				subvector(ub_, size_.nu(), size_.nx()) = ubx; 
			}

			auto const& x() const { return x_; }
			auto const& u() const { return u_;	}
			auto const& pi() const	{ return pi_; }
			
			auto lam_lbu() const 
			{ 
				return subvector(lam_, 2 * size_.nc(), size_.nu()); 
			}

			auto lam_ubu() const 
			{ 
				return subvector(lam_, 2 * size_.nc() + size_.nu(), size_.nu()); 
			}

			auto lam_lbx() const 
			{ 
				return subvector(lam_, 2 * size_.nc() + 2 * size_.nu(), size_.nx()); 
			}

			auto lam_ubx() const 
			{ 
				return subvector(lam_, 2 * size_.nc() + 2 * size_.nu() + size_.nx(), size_.nx()); 
			}

			auto lam_lbd() const 
			{ 
				return subvector(lam_, 0, size_.nc()); 
			}

			auto lam_ubd() const 
			{ 
				return subvector(lam_, size_.nc(), size_.nc()); 
			}

			QpSize const& size() const { return size_; }

			// Adjust hidxb so to account for infs in state and input bounds.
			void adjustBoundsIndex()
			{
				// this will not change the capacity and the data() pointers should stay the same.
				hidxb_.clear();
				lb_internal_.clear();
				ub_internal_.clear();

				// Cycle through the bounds and check for infinities
				for (size_t i = 0; i < size_.nu() + size_.nx(); ++i)
				{
					if (std::isfinite(lb_[i]) && std::isfinite(ub_[i]))
					{
						// If both bounds are finite, add i to the bounds index,
						// and copy values to the lb_internal_ and ub_internal_.
						hidxb_.push_back(i);
						lb_internal_.push_back(lb_[i]);
						ub_internal_.push_back(ub_[i]);
					}
					else 
					{
						// Otherwise, check that the values are [-inf, inf]
						if (!(lb_[i] == -inf<Real>() && ub_[i] == inf<Real>()))
							throw std::invalid_argument("And invalid QP bound is found. For HPMPC, "
								"the bounds should be either both finite or [-inf, inf]");
					}
				}
			}

			// ******************************************************
			//                HPMPC raw data interface.
			//
			// The prefixes before _data() correspond to the names of
			// the argument to c_order_d_ip_ocp_hard_tv().
			// ******************************************************
			Real const * A_data () const { return A_.data(); }
			Real const * B_data () const { return B_.data(); }
			Real const * b_data () const { return b_.data();	}
			Real const * Q_data () const { return Q_.data(); }
			Real const * S_data () const { return S_.data(); }
			Real const * R_data () const { return R_.data(); }
			Real const * q_data () const { return q_.data(); }
			Real const * r_data () const { return r_.data();	}
			Real const * lb_data() const { return lb_internal_.data(); }
			Real const * ub_data() const { return ub_internal_.data(); }
			Real const * C_data () const { return C_.data(); }
			Real const * D_data () const { return D_.data(); }
			Real const * lg_data() const { return lbd_.data(); }
			Real const * ug_data() const { return ubd_.data(); }
			int const * hidxb_data() const { return hidxb_.data(); }
			int nb() const { return hidxb_.size(); }

			Real * x_data() { return x_.data(); }
			Real * u_data() { return u_.data(); }
			Real * pi_data() { return pi_.data(); }
			Real * lam_data() { return lam_.data(); }

		private:
			QpSize size_;

			// Some magic data for HPMPC
			std::vector<int> hidxb_;

			// Hessian = [R, S; S', Q]
			Matrix R_;
			Matrix S_;
			Matrix Q_;

			// Gradient = [r; q]
			Vector r_;
			Vector q_;

			// Inter-stage equalities x_{k+1} = A x_k + B u_k + c_k
			Matrix A_;
			Matrix B_;
			Vector b_;

			// Inequality constraints d_{min} <= C x_k + D u_k <= d_{max}
			Matrix C_;
			Matrix D_;
			Vector lbd_;
			Vector ubd_;

			// Bound constraints:
			// lb <= [u; x] <= ub
			Vector lb_;
			Vector ub_;

			// Lower and upper bound arrays for HPMPC,
			// containing finite values only.
			std::vector<Real> lb_internal_;
			std::vector<Real> ub_internal_;

			Vector x_;
			Vector u_;
			Vector pi_;
			Vector lam_;
		};

		class ProblemIterator
		:	public boost::iterator_adaptor<
				ProblemIterator	// derived
			,	typename std::vector<Stage>::iterator	// base
			,	QpStageBase<Stage>&	// value
			,	boost::random_access_traversal_tag	// category of traversal
			,	QpStageBase<Stage>&	// reference
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
			,	QpStageBase<Stage> const&	// value
			,	boost::random_access_traversal_tag	// category of traversal
			,	QpStageBase<Stage> const&	// reference
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
			,	QpStageSolutionBase<Stage> const&	// value
			,	boost::random_access_traversal_tag	// category of traversal
			,	QpStageSolutionBase<Stage> const&	// reference
			>
		{
		public:
			ConstSolutionIterator() = default;
			
			ConstSolutionIterator(typename ConstSolutionIterator::base_type const& p)
			:	ConstSolutionIterator::iterator_adaptor_(p)
			{			
			}
		};
	
		boost::iterator_range<ProblemIterator> problem()
		{
			return boost::iterator_range<ProblemIterator>(stage_.begin(), stage_.end());
		}
	
		boost::iterator_range<ConstProblemIterator> problem() const
		{
			return boost::iterator_range<ConstProblemIterator>(stage_.begin(), stage_.end());
		}
	
		boost::iterator_range<ConstSolutionIterator> solution() const
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
		explicit HpmpcWorkspace(InputIterator size_first, InputIterator size_last, int max_iter = 100)
		:	stat_(max_iter)
		,	infNormRes_(sNaN<Real>())
		{
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

		explicit HpmpcWorkspace(std::initializer_list<QpSize> sz, int max_iter = 100)
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

		void solve()
		{
			if (stage_.size() > 0)
			{
				// Recalculate bounds indices of each stage, at the bounds might have changed.				
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
		typename Kernel::template StaticVector<4> infNormRes_;

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
		void addStage(QpSize const& sz, size_t nx_next)
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
