#pragma once

#include "QpSolverException.hpp"
#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpSolutionBase.hpp>
#include <tmpc/qp/OcpQpBase.hpp>
#include <tmpc/qp/QpWorkspaceBase.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

#include "detail/HpxxxStage.hpp"

#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_ipm.h>

#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>


namespace tmpc
{
	template <typename Real>
	struct Hpipm;

	template <>
	struct Hpipm<double>
	{
		using ocp_qp = ::d_ocp_qp;
		using ocp_qp_sol = ::d_ocp_qp_sol;
		using ocp_qp_ipm_workspace = ::d_ocp_qp_ipm_workspace;
		using ocp_qp_ipm_arg = ::d_ocp_qp_ipm_arg;

		static int memsize_ocp_qp(int N, int const *nx, int const *nu, int const *nb, int const *ng, int const * ns);
		static void create_ocp_qp(int N, int const *nx, int const *nu, int const *nb, int const *ng, int const * ns, ocp_qp *qp, void *memory);
		static void cvt_colmaj_to_ocp_qp(
			double const * const *A, double const * const *B, double const * const *b, 
			double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r, 
			int const * const *idxb, double const * const *lb, double const * const *ub, 
			double const * const *C, double const * const *D, double const * const *lg, double const * const *ug,
			double const * const *Zl, double const * const *Zu, double const * const *zl, double const * const *zu, int const * const *idxs, 
			ocp_qp * qp);
		
		static int memsize_ocp_qp_sol(int N, int const *nx, int const *nu, int const *nb, int const *ng, int const *ns);
		static void create_ocp_qp_sol(int N, int const *nx, int const *nu, int const *nb, int const *ng, int const *ns, ocp_qp_sol *qp_sol, void *memory);
		static void cvt_ocp_qp_sol_to_colmaj(ocp_qp const *qp, ocp_qp_sol const *qp_sol, 
			double * const *u, double * const *x, double * const * ls, double * const * us,
			double * const *pi, 
			double * const *lam_lb, double * const *lam_ub, double * const *lam_lg, double * const *lam_ug,
			double * const *lam_ls, double * const *lam_us);

		static int memsize_ocp_qp_ipm(ocp_qp const *qp, ocp_qp_ipm_arg const *arg);
		static void create_ocp_qp_ipm(ocp_qp const *qp, ocp_qp_ipm_arg const *arg, ocp_qp_ipm_workspace *ws, void *mem);
		static int solve_ocp_qp_ipm(ocp_qp const *qp, ocp_qp_sol *qp_sol, ocp_qp_ipm_arg const *arg, ocp_qp_ipm_workspace *ws);
	};

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
	:	public QpWorkspaceBase<HpipmWorkspace<Kernel_>>
	{
	public:
		using Kernel = Kernel_;
		using Real = typename Kernel::Real;
		
	private:
		using Stage = detail::HpxxxStage<Kernel, StorageOrder::columnMajor>;

	public:

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
			return "HPIPM";
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

		/**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <typename InputIterator>
		HpipmWorkspace(InputIterator size_first, InputIterator size_last, int max_iter = 100)
		:	solverArg_ {}
		{
			std::fill(infNormRes_.begin(), infNormRes_.end(), sNaN());

			solverArg_.alpha_min = 1e-8;
			solverArg_.res_g_max = 1e-8;
			solverArg_.res_b_max = 1e-8;
			solverArg_.res_d_max = 1e-12;
			solverArg_.res_m_max = 1e-12;
			solverArg_.mu0 = 2.0;
			solverArg_.iter_max = max_iter;
			solverArg_.stat_max = 100;
			solverArg_.pred_corr = 1;

			auto const nt = std::distance(size_first, size_last);
			
			// Preallocate arrays holding QP stage data.
			stage_.reserve(nt);
	
			nx_.reserve(nt);
			nu_.reserve(nt);
			nb_.reserve(nt);
			ng_.reserve(nt);
			ns_.reserve(nt);
			hidxb_.reserve(nt);
			idxs_.reserve(nt);
	
			A_ .reserve(nt);
			B_ .reserve(nt);
			b_ .reserve(nt);
			Q_ .reserve(nt);
			S_ .reserve(nt);
			R_ .reserve(nt);
			q_ .reserve(nt);
			r_ .reserve(nt);
			Zl_.reserve(nt);
			Zu_.reserve(nt);
			zl_.reserve(nt);
			zu_.reserve(nt);
			lb_.reserve(nt);
			ub_.reserve(nt);
			C_ .reserve(nt);
			D_ .reserve(nt);
			lg_.reserve(nt);
			ug_.reserve(nt);
			
			x_.reserve(nt);
			u_.reserve(nt);
			pi_.reserve(nt);
			lam_lb_.reserve(nt);
			lam_ub_.reserve(nt);
			lam_lg_.reserve(nt);
			lam_ug_.reserve(nt);

			for (auto sz = size_first; sz != size_last; ++sz)
			{
				auto const nx_next = sz + 1 != size_last ? sz[1].nx() : 0;

				stage_.emplace_back(*sz, nx_next);
				auto& st = stage_.back();
		
				nx_.push_back(sz->nx());
				nu_.push_back(sz->nu());
				nb_.push_back(st.nb());
				ng_.push_back(sz->nc());
				ns_.push_back(sz->ns());
				
				hidxb_.push_back(st.hidxb_data());
				idxs_.push_back(st.idxs_data());
		
				A_.push_back(st.A_data());
				B_.push_back(st.B_data());
				b_.push_back(st.b_data());
		
				Q_.push_back(st.Q_data());
				S_.push_back(st.S_data());
				R_.push_back(st.R_data());
				q_.push_back(st.q_data());
				r_.push_back(st.r_data());
				Zl_.push_back(st.Zl_data());
				Zu_.push_back(st.Zu_data());
				zl_.push_back(st.zl_data());
				zu_.push_back(st.zu_data());
				lb_.push_back(st.lb_data());
				ub_.push_back(st.ub_data());
		
				C_ .push_back(st.C_data());
				D_ .push_back(st.D_data());
				lg_.push_back(st.lg_data());
				ug_.push_back(st.ug_data());
		
				x_.push_back(st.x_data());
				u_.push_back(st.u_data());
				ls_.push_back(st.ls_data());
				us_.push_back(st.us_data());
				pi_.push_back(st.pi_data());
				lam_lb_.push_back(st.lam_lb_data());
				lam_ub_.push_back(st.lam_ub_data());
				lam_lg_.push_back(st.lam_lg_data());
				lam_ug_.push_back(st.lam_ug_data());
				lam_ls_.push_back(st.lam_ls_data());
				lam_us_.push_back(st.lam_us_data());
			}
				
			if (stage_.size() > 1)
			{
				int const N = stage_.size() - 1;

				// Allocate big enough memory pools for Qp, QpSol and SolverWorkspace.
				// nb_ is set to nx+nu at this point, which ensures maximum capacity.
				ocpQpMem_.resize(HPIPM::memsize_ocp_qp(N, nx_.data(), nu_.data(), nb_.data(), ng_.data(), ns_.data()));

				// This ocp_qp borns and dies just for memsize_ipm_hard_ocp_qp() to calculate the necessary workspace size.
				typename HPIPM::ocp_qp qp;
				HPIPM::create_ocp_qp(N, nx_.data(), nu_.data(), nb_.data(), ng_.data(), ns_.data(), &qp, ocpQpMem_.data());

				ocpQpSolMem_.resize(HPIPM::memsize_ocp_qp_sol(N, nx_.data(), nu_.data(), nb_.data(), ng_.data(), ns_.data()));
				solverWorkspaceMem_.resize(HPIPM::memsize_ocp_qp_ipm(&qp, &solverArg_));				
			}
		}

		template <typename IteratorRange>
		explicit HpipmWorkspace(IteratorRange sz, int max_iter = 100)
		:	HpipmWorkspace(sz.begin(), sz.end(), max_iter)
		{
		}

		explicit HpipmWorkspace(std::initializer_list<OcpSize> sz, int max_iter = 100)
		:	HpipmWorkspace(sz.begin(), sz.end(), max_iter)
		{
		}

		explicit HpipmWorkspace(std::initializer_list<OcpSize> sz, int max_iter = 100)
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

		void impl_solve()
		{
			if (stage_.size() > 1)
			{
				// Recalculate bounds indices of each stage, at the bounds might have changed.				
				// Update the nb_ array.
				std::transform(stage_.begin(), stage_.end(), nb_.begin(), [] (Stage& s) -> int
				{
					s.adjustBoundsIndex();
					return s.nb();
				});
	
				// Number of QP steps for HPIPM
				auto const N = stage_.size() - 1;
	
				// Init HPIPM structures. Since the nb_ might change if the bounds change, we need to do it every time, sorry.
				// Hopefully these operations are not expensive compared to actual solving.
				typename HPIPM::ocp_qp qp;
				HPIPM::create_ocp_qp(N, nx_.data(), nu_.data(), nb_.data(), ng_.data(), ns_.data(), &qp, ocpQpMem_.data());
	
				typename HPIPM::ocp_qp_sol sol;
				HPIPM::create_ocp_qp_sol(N, nx_.data(), nu_.data(), nb_.data(), ng_.data(), ns_.data(), &sol, ocpQpSolMem_.data());
	
				NativeHpipmWorkspace solver_workspace;
				HPIPM::create_ocp_qp_ipm(&qp, &solverArg_, &solver_workspace, solverWorkspaceMem_.data());
	
				// Convert the problem
				HPIPM::cvt_colmaj_to_ocp_qp(
					A_.data(), B_.data(), b_.data(), 
					Q_.data(), S_.data(), R_.data(), q_.data(), r_.data(), 
					hidxb_.data(), lb_.data(), ub_.data(), 
					C_.data(), D_.data(), lg_.data(), ug_.data(), 
					Zl_.data(), Zu_.data(), zl_.data(), zu_.data(), idxs_.data(),
					&qp);
	
				// Call HPIPM
				auto const ret = HPIPM::solve_ocp_qp_ipm(&qp, &sol, &solverArg_, &solver_workspace);
				numIter_ = solver_workspace.iter;
	
				if (ret != 0)
				{
					throw HpipmException(ret);
				}
	
				// Convert the solution
				HPIPM::cvt_ocp_qp_sol_to_colmaj(&qp, &sol, 
					u_.data(), x_.data(), ls_.data(), us_.data(),
					pi_.data(), lam_lb_.data(), lam_ub_.data(), lam_lg_.data(), lam_ug_.data(),
					lam_ls_.data(), lam_us_.data());
			}
		}

		std::size_t maxIter() const noexcept { return solverArg_.iter_max; }

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

		
	private:
		using HPIPM = Hpipm<Real>;
		using NativeHpipmWorkspace = typename HPIPM::ocp_qp_ipm_workspace;

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

		std::vector<Real const *> Zl_;
		std::vector<Real const *> Zu_;
		std::vector<Real const *> zl_;
		std::vector<Real const *> zu_;

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

		// Array of NS (soft constraints) sizes
		std::vector<int> ns_;

		// Hard constraints index
		std::vector<int const *> hidxb_;

		// Soft constraints index
		std::vector<int const *> idxs_;

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

		typename HPIPM::ocp_qp_ipm_arg solverArg_;


		// Warmstarting disabled on purpose.
		// On AMD K8 (hpmpc compiled for SSE3), WITHOUT warmstarting it is significantly
		// FASTER (9ms vs 14ms per time step) than with warmstarting. I am curious why.
		bool _warmStart = false;

		static Real constexpr sNaN()
		{
			return std::numeric_limits<Real>::signaling_NaN();
		}
	};
}
