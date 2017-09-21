#pragma once

#include "QpSolverException.hpp"
#include "QpSize.hpp"
#include "QpStageSolutionBase.hpp"
#include "QpStageBase.hpp"

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_ipm_hard.h>

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>

namespace tmpc
{
	template <typename Real>
	struct Hpipm;

	template <>
	struct Hpipm<double>
	{
		using ocp_qp = ::d_ocp_qp;
		using ocp_qp_sol = ::d_ocp_qp_sol;
		using ipm_hard_ocp_qp_workspace = ::d_ipm_hard_ocp_qp_workspace;
		using ipm_hard_ocp_qp_arg = ::d_ipm_hard_ocp_qp_arg;

		static int memsize_ocp_qp(int N, int const *nx, int const *nu, int const *nb, int const *ng);
		static void create_ocp_qp(int N, int const *nx, int const *nu, int const *nb, int const *ng, ocp_qp *qp, void *memory);
		static void cvt_colmaj_to_ocp_qp(
			double const * const *A, double const * const *B, double const * const *b, 
			double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r, 
			int const * const *idxb, double const * const *lb, double const * const *ub, 
			double const * const *C, double const * const *D, double const * const *lg, double const * const *ug, ocp_qp * qp);
		
		static int memsize_ocp_qp_sol(int N, int const *nx, int const *nu, int const *nb, int const *ng);
		static void create_ocp_qp_sol(int N, int const *nx, int const *nu, int const *nb, int const *ng, ocp_qp_sol *qp_sol, void *memory);
		static void cvt_ocp_qp_sol_to_colmaj(ocp_qp const *qp, ocp_qp_sol const *qp_sol, 
			double * const *u, double * const *x, double * const *pi, double * const *lam_lb, double * const *lam_ub, double * const *lam_lg, double * const *lam_ug);

		static int memsize_ipm_hard_ocp_qp(ocp_qp const *qp, ipm_hard_ocp_qp_arg const *arg);
		static void create_ipm_hard_ocp_qp(ocp_qp const *qp, ipm_hard_ocp_qp_arg const *arg, ipm_hard_ocp_qp_workspace *ws, void *mem);
		static void solve_ipm_hard_ocp_qp(ocp_qp const *qp, ocp_qp_sol *qp_sol, ipm_hard_ocp_qp_workspace *ws);
	};

	class HpipmException 
	:	public QpSolverException
	{
	public:
		HpipmException(int code);

		int code() const { return _code; }
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
	class HpipmWorkspace
	{
	public:
		using Real = Real_;

		class Stage
		:	public QpStageSolutionBase<Stage>
		,	public QpStageBase<Stage>
		{
		public:
			static auto constexpr storageOrder = columnMajor;
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
			template <typename T> void A(const T& a) { noresize(A_) = a; }

			const Matrix& B() const { return B_; }
			template <typename T> void B(const T& b) { noresize(B_) = b; }

			Vector const& b() const { return b_; }
			template <typename T> void b(const T& b) { noresize(b_) = b; }

			const Matrix& C() const { return C_; }
			template <typename T> void C(const T& c) { noresize(C_) = c; }

			const Matrix& D() const { return D_; }
			template <typename T> void D(const T& d) { noresize(D_) = d; }

			const Vector& lbd() const {	return lbd_; }
			template <typename T> void lbd(const T& lbd) { noresize(lbd_) = lbd; }

			auto lbu() const { return subvector(lb_, 0, size_.nu());	}			
			template <typename T> void lbu(const T& lbu) { subvector(lb_, 0, size_.nu()) = lbu; }

			auto lbx() const { return subvector(lb_, size_.nu(), size_.nx()); }
			template <typename T> void lbx(const T& lbx) { subvector(lb_, size_.nu(), size_.nx()) = lbx; }

			const Matrix& Q() const { return Q_; }
			template <typename T> void Q(const T& q) { noresize(Q_) = q; }

			const Matrix& R() const { return R_; }
			template <typename T> void R(const T& r) { noresize(R_) = r; }

			// HPIPM convention for S is [nu, nx], therefore the trans().
			decltype(auto) S() const { return trans(S_); }
			template <typename T> void S(const T& s) { noresize(S_) = trans(s); }

			const Vector& q() const { return q_; }
			template <typename T> void q(const T& q) { noresize(q_) = q; }

			const Vector& r() const { return r_; }
			template <typename T> void r(const T& r) { noresize(r_) = r; }

			const Vector& ubd() const { return ubd_; }
			template <typename T> void ubd(const T& ubd) { noresize(ubd_) = ubd; }

			auto ubu() const { return subvector(ub_, 0, size_.nu()); }
			template <typename T> void ubu(const T& ubu) { subvector(ub_, 0, size_.nu()) = ubu; }

			auto ubx() const { return subvector(ub_, size_.nu(), size_.nx()); }
			template <typename T> void ubx(const T& ubx) { subvector(ub_, size_.nu(), size_.nx()) = ubx; }

			DynamicVector<Real> const& x() const { return x_; }
			DynamicVector<Real> const& u() const { return u_;	}
			DynamicVector<Real> const& pi() const	{ return pi_; }
			
			decltype(auto) lam_lbu() const 
			{ 
				return subvector(lam_, 0, size_.nu()); 
			}

			decltype(auto) lam_ubu() const 
			{ 
				return subvector(lam_, size_.nu(), size_.nu()); 
			}

			decltype(auto) lam_lbx() const 
			{ 
				return subvector(lam_, 2 * size_.nu(), size_.nx()); 
			}

			decltype(auto) lam_ubx() const 
			{ 
				return subvector(lam_, 2 * size_.nu() + size_.nx(), size_.nx()); 
			}

			decltype(auto) lam_lbd() const 
			{ 
				return subvector(lam_, 2 * size_.nu() + 2 * size_.nx(), size_.nc()); 
			}

			decltype(auto) lam_ubd() const 
			{ 
				return subvector(lam_, 2 * size_.nu() + 2 * size_.nx() + size_.nc(), size_.nc()); 
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
			//                HPIPM raw data interface.
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
			Real * lam_lb_data() { return lam_.data(); }
			Real * lam_ub_data() { return lam_lb_data() + size_.nx() + size_.nu(); }
			Real * lam_lg_data() { return lam_ub_data() + size_.nx() + size_.nu(); }
			Real * lam_ug_data() { return lam_lg_data() + size_.nc(); }

		private:
			QpSize size_;

			// Some magic data for HPIPM
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

			// Lower and upper bound arrays for HPIPM,
			// containing finite values only.
			std::vector<Real> lb_internal_;
			std::vector<Real> ub_internal_;

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

		/**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <typename InputIterator>
		HpipmWorkspace(InputIterator size_first, InputIterator size_last, int max_iter = 100)
		:	infNormRes_(sNaN())
		,	solverArg_ {}
		{
			solverArg_.alpha_min = 1e-8;
			solverArg_.mu_max = 1e-12;
			solverArg_.iter_max = max_iter;
			solverArg_.mu0 = 2.0;
			
			preallocateStages(std::distance(size_first, size_last));

			for (auto sz = size_first; sz != size_last; ++sz)
				addStage(*sz, sz + 1 != size_last ? sz[1].nx() : 0);
				
			if (stage_.size() > 1)
			{
				int const N = stage_.size() - 1;

				// Allocate big enough memory pools for Qp, QpSol and SolverWorkspace.
				// nb_ is set to nx+nu at this point, which ensures maximum capacity.
				ocpQpMem_.resize(HPIPM::memsize_ocp_qp(N, nx_.data(), nu_.data(), nb_.data(), ng_.data()));

				// This ocp_qp borns and dies just for memsize_ipm_hard_ocp_qp() to calculate the necessary workspace size.
				typename HPIPM::ocp_qp qp;
				HPIPM::create_ocp_qp(N, nx_.data(), nu_.data(), nb_.data(), ng_.data(), &qp, ocpQpMem_.data());

				ocpQpSolMem_.resize(HPIPM::memsize_ocp_qp_sol(N, nx_.data(), nu_.data(), nb_.data(), ng_.data()));
				solverWorkspaceMem_.resize(HPIPM::memsize_ipm_hard_ocp_qp(&qp, &solverArg_));				
			}
		}

		template <typename IteratorRange>
		explicit HpipmWorkspace(IteratorRange sz, int max_iter = 100)
		:	HpipmWorkspace(sz.begin(), sz.end(), max_iter)
		{
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

		void solve();

		std::size_t maxIter() const noexcept { return solverArg_.iter_max; }

		/// \brief Get number of iterations performed by the QP solver.
		unsigned numIter() const { return numIter_; }
		
	private:
		using HPIPM = Hpipm<Real>;

		class OcpQp
		{
		public:
			OcpQp();
			OcpQp(int N, int const *nx, int const *nu, int const *nb, int const *ng);
	
		private:
			std::vector<char> memory_;
			typename HPIPM::ocp_qp qp_;
		};

		class OcpQpSol
		{
		public:
			OcpQpSol(int N, int const *nx, int const *nu, int const *nb, int const *ng);
	
		private:
			std::vector<char> memory_;
			typename HPIPM::ocp_qp_sol sol_;
		};

		class IpmHardOcpQpWorkspace
		{
		public:
			IpmHardOcpQpWorkspace(typename HPIPM::ocp_qp const& qp, typename HPIPM::ipm_hard_ocp_qp_arg const& arg);
	
		private:
			typename HPIPM::ipm_hard_ocp_qp_workspace workspace_;
			std::vector<char> memory_;
		};

		// --------------------------------
		//
		// QP problem data
		//
		// --------------------------------

		// Stores stage data
		std::vector<Stage> stage_;

		// "A" data array for HPIPM
		std::vector<Real const *> A_;

		// "B" data array for HPIPM
		std::vector<Real const *> B_;

		// "b" data array for HPIPM
		std::vector<Real const *> b_;

		// "Q" data array for HPIPM
		std::vector<Real const *> Q_;

		// "S" data array for HPIPM
		std::vector<Real const *> S_;

		// "R" data array for HPIPM
		std::vector<Real const *> R_;

		// "q" data array for HPIPM
		std::vector<Real const *> q_;

		// "r" data array for HPIPM
		std::vector<Real const *> r_;

		// "lb" data array for HPIPM
		std::vector<Real const *> lb_;

		// "ub" data array for HPIPM
		std::vector<Real const *> ub_;

		// "C" data array for HPIPM
		std::vector<Real const *> C_;

		// "D" data array for HPIPM
		std::vector<Real const *> D_;

		// "lg" data array for HPIPM
		std::vector<Real const *> lg_;

		// "ug" data array for HPIPM
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
		// QP solution data
		//
		// --------------------------------

		std::vector<Real *> x_;
		std::vector<Real *> u_;
		std::vector<Real *> pi_;
		std::vector<Real *> lam_lb_;
		std::vector<Real *> lam_ub_;
		std::vector<Real *> lam_lg_;
		std::vector<Real *> lam_ug_;
		StaticVector<Real, 4> infNormRes_;

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

		typename HPIPM::ipm_hard_ocp_qp_arg solverArg_;


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
	};
}
