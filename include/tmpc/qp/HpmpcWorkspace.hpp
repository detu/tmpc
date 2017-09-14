#pragma once

#include "QpSolverException.hpp"
#include "QpSize.hpp"
#include "QpStageSolutionBase.hpp"
#include "QpStageBase.hpp"

#include <tmpc/Matrix.hpp>

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>

namespace tmpc
{
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
	 * \brief Multistage QP solver using qpOASES
	 *
	 * \tparam <Real_> real number type
	 */
	template <typename Real_>
	class HpmpcWorkspace
	{
	public:
		using Real = Real_;
		
		// Iteration statistics. HPMPC returns 5 Real numbers per iteration.
		using IterStat = std::array<Real, 5>;

		class Stage
		:	public QpStageSolutionBase<Stage>
		,	public QpStageBase<Stage>
		{
		public:
			static auto constexpr storageOrder = rowMajor;
			using Real = HpmpcWorkspace::Real;
			typedef DynamicMatrix<Real, storageOrder> Matrix;
			typedef DynamicVector<Real, columnVector> Vector;

			Stage(QpSize const& sz, size_t nx_next);
			Stage(Stage const&) = default;
			Stage(Stage &&) = default;

			template <typename Expr>
			Stage& operator=(Expr const& rhs)
			{
				assign(*this, rhs);
				return *this;
			}

			const Matrix& A() const { return A_; }
			template <typename T> void A(const T& a) { full(A_) = a; }

			const Matrix& B() const { return B_; }
			template <typename T> void B(const T& b) { full(B_) = b; }

			Vector const& b() const { return b_; }
			template <typename T> void b(const T& b) { full(b_) = b; }

			const Matrix& C() const { return C_; }
			template <typename T> void C(const T& c) { full(C_) = c; }

			const Matrix& D() const { return D_; }
			template <typename T> void D(const T& d) { full(D_) = d; }

			const Vector& lbd() const {	return lbd_; }
			template <typename T> void lbd(const T& lbd) { full(lbd_) = lbd; }

			Subvector<Vector const> lbu() const { return subvector(lb_, 0, size_.nu());	}			
			template <typename T> void lbu(const T& lbu) { subvector(lb_, 0, size_.nu()) = lbu; }

			Subvector<Vector const> lbx() const { return subvector(lb_, size_.nu(), size_.nx()); }
			template <typename T> void lbx(const T& lbx) { subvector(lb_, size_.nu(), size_.nx()) = lbx; }

			const Matrix& Q() const { return Q_; }
			template <typename T> void Q(const T& q) { full(Q_) = q; }

			const Matrix& R() const { return R_; }
			template <typename T> void R(const T& r) { full(R_) = r; }

			// HPMPC convention for S is [nu, nx], therefore the trans().
			decltype(auto) S() const { return trans(S_); }
			template <typename T> void S(const T& s) { full(S_) = trans(s); }

			const Vector& q() const { return q_; }
			template <typename T> void q(const T& q) { full(q_) = q; }

			const Vector& r() const { return r_; }
			template <typename T> void r(const T& r) { full(r_) = r; }

			const Vector& ubd() const { return ubd_; }
			template <typename T> void ubd(const T& ubd) { full(ubd_) = ubd; }

			Subvector<Vector const> ubu() const { return subvector(ub_, 0, size_.nu()); }
			template <typename T> void ubu(const T& ubu) { subvector(ub_, 0, size_.nu()) = ubu; }

			Subvector<Vector const> ubx() const { return subvector(ub_, size_.nu(), size_.nx()); }
			template <typename T> void ubx(const T& ubx) { subvector(ub_, size_.nu(), size_.nx()) = ubx; }

			DynamicVector<Real> const& x() const { return x_; }
			DynamicVector<Real> const& u() const { return u_;	}
			DynamicVector<Real> const& pi() const	{ return pi_; }
			
			decltype(auto) lam_lbu() const 
			{ 
				return subvector(lam_, 2 * size_.nc(), size_.nu()); 
			}

			decltype(auto) lam_ubu() const 
			{ 
				return subvector(lam_, 2 * size_.nc() + size_.nu(), size_.nu()); 
			}

			decltype(auto) lam_lbx() const 
			{ 
				return subvector(lam_, 2 * size_.nc() + 2 * size_.nu(), size_.nx()); 
			}

			decltype(auto) lam_ubx() const 
			{ 
				return subvector(lam_, 2 * size_.nc() + 2 * size_.nu() + size_.nx(), size_.nx()); 
			}

			decltype(auto) lam_lbd() const 
			{ 
				return subvector(lam_, 0, size_.nc()); 
			}

			decltype(auto) lam_ubd() const 
			{ 
				return subvector(lam_, size_.nc(), size_.nc()); 
			}

			QpSize const& size() const { return size_; }

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
			Real const * lb_data() const { return lb_.data(); }
			Real const * ub_data() const { return ub_.data(); }
			Real const * C_data () const { return C_.data(); }
			Real const * D_data () const { return D_.data(); }
			Real const * lg_data() const { return lbd_.data(); }
			Real const * ug_data() const { return ubd_.data(); }
			int const * hidxb_data() const { return hidxb_.data(); }

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

			DynamicVector<Real> x_;
			DynamicVector<Real> u_;
			DynamicVector<Real> pi_;
			DynamicVector<Real> lam_;
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
		HpmpcWorkspace(InputIterator size_first, InputIterator size_last, int max_iter = 100)
		:	stat_(max_iter)
		,	infNormRes_(sNaN())
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

		void solve();

		std::size_t maxIter() const noexcept { return stat_.size(); }

		/// \brief Set max number of iterations
		void maxIter(size_t val);

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
		// --------------------------------
		//
		// HPMPC QP problem data
		//
		// --------------------------------

		// Stores stage data
		std::vector<Stage> stage_;

		// "A" data array for HPMPC
		std::vector<Real const *> A_;

		// "B" data array for HPMPC
		std::vector<Real const *> B_;

		// "b" data array for HPMPC
		std::vector<Real const *> b_;

		// "Q" data array for HPMPC
		std::vector<Real const *> Q_;

		// "S" data array for HPMPC
		std::vector<Real const *> S_;

		// "R" data array for HPMPC
		std::vector<Real const *> R_;

		// "q" data array for HPMPC
		std::vector<Real const *> q_;

		// "r" data array for HPMPC
		std::vector<Real const *> r_;

		// "lb" data array for HPMPC
		std::vector<Real const *> lb_;

		// "ub" data array for HPMPC
		std::vector<Real const *> ub_;

		// "C" data array for HPMPC
		std::vector<Real const *> C_;

		// "D" data array for HPMPC
		std::vector<Real const *> D_;

		// "lg" data array for HPMPC
		std::vector<Real const *> lg_;

		// "ug" data array for HPMPC
		std::vector<Real const *> ug_;

		// Array of NX sizes
		std::vector<int> nx_;

		// Array of NU sizes
		std::vector<int> nu_;

		// Array of NB (bound constraints) sizes
		std::vector<int> nb_;

		// Array of NG (path constraints) sizes
		std::vector<int> ng_;

		// Additional data for new hpmpc interface
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
		StaticVector<Real, 4> infNormRes_;

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

		static Real constexpr sNaN()
		{
			return std::numeric_limits<Real>::signaling_NaN();
		}

		// Preallocate arrays holding QP stage data.
		void preallocateStages(size_t nt);

		// Add one stage with specified sizes at the end of the stage sequence.
		void addStage(QpSize const& sz, size_t nx_next);

		// Allocate soverWorkspace_ according to nx, nu, nb, ng etc.
		void allocateSolverWorkspace();
	};
}
