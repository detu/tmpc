#pragma once

#include "HPMPCProblem.hpp"
#include "HPMPCSolution.hpp"
#include "UnsolvedQpException.hpp"

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>


namespace tmpc
{
	namespace hpmpc_wrapper
	{	
		int c_order_d_ip_ocp_hard_tv(	int *kk, int k_max, double mu0, double mu_tol,
										int N, int const *nx, int const *nu, int const *nb, int **hidxb, int const *ng, int N2, 
										int warm_start,
										double const * const *A, double const * const *B, double const * const *b,
										double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r,
										double const * const *lb, double const * const *ub,
										double const * const *C, double const * const *D, double const * const *lg, double const * const *ug,
										double * const *x, double * const *u, double * const *pi, double * const *lam, double * const *t,
										double *inf_norm_res,
										void *work0,
										double *stat );

		int d_ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int **hidxb, int const *ng, int N2);
		
	}

	class HpmpcUnsolvedQpException : public UnsolvedQpException
	{
	public:
		template <typename QP>
		HpmpcUnsolvedQpException(QP const& qp, int code)
		:	UnsolvedQpException("HPMPC", qp),
			_code(code),
			msg_(std::string(UnsolvedQpException::what()) + "\nReturn code " + std::to_string(code))
		{
		}

		int getCode() const	{ return _code;	}
		char const * what() const noexcept override { return msg_.c_str(); }

	private:
		int const _code;
		std::string const msg_;
	};

	/**
	 * \brief Multistage QP solver using qpOASES
	 *
	 * \tparam <K> Class implementing the Kernel concept
	 * \tparam <D> Class defining problem dimensions
	 */
	template <typename K, typename D>
	class HPMPCSolver
	{
		static auto constexpr NX = D::NX;
		static auto constexpr NU = D::NU;
		static auto constexpr NZ = D::NX + D::NU;
		static auto constexpr NC = D::NC;
		static auto constexpr NCT = D::NCT;

		typedef typename K::Scalar Scalar;

	public:

		typedef HPMPCProblem<NX, NU, NC, NCT> Problem;
		typedef HPMPCSolution<NX, NU, NC, NCT> Solution;

		HPMPCSolver(std::size_t nt, int max_iter = 100)
		:	_nx(nt + 1, NX)
		,	_nu(nt + 1, NU)
		,	_nb(nt + 1, NU + NX)
		,	_ng(nt + 1, NC)
		,	_stat(max_iter)
		,	stageBounds_(nt)
		,	lb_(nt + 1)
		,	ub_(nt + 1)
		,	lg_(nt + 1)
		,	ug_(nt + 1)
		{
			_nu.back() = 0;
			_nb.back() = NX;
			_ng.back() = NCT;

			for (std::size_t i = 0; i < nt; ++i)
			{
				lb_[i] = stageBounds_[i].lb.data();
				ub_[i] = stageBounds_[i].ub.data();
				lg_[i] = stageBounds_[i].lg.data();
				ug_[i] = stageBounds_[i].ug.data();
			}

			lb_.back() = terminalStageBounds_.lb.data();
			ub_.back() = terminalStageBounds_.ub.data();
			lg_.back() = terminalStageBounds_.lg.data();
			ug_.back() = terminalStageBounds_.ug.data();

			// Data structure for new hpmpc interface
			_hidxb = new int *[nt+1];
			for(int kk = 0; kk < nt; kk++) {
				_hidxb[kk] = new int[NU + NX];
				for (int ii = 0; ii < NU + NX; ii++) {
					_hidxb[kk][ii] = ii;	
				}	
			}
			
			_hidxb[nt] = new int[NX];
			for (int ii = 0; ii < NX; ii++) {
				_hidxb[nt][ii] = ii;
			}

			// Allocate workspace
			_workspace.resize(hpmpc_wrapper::d_ip_ocp_hard_tv_work_space_size_bytes(
					static_cast<int>(nt), _nx.data(), _nu.data(), _nb.data(), _hidxb, _ng.data(), static_cast<int>(nt)));
		}
		
		/**
		 * \brief Destructor
		 */
		~HPMPCSolver() {
			if(_hidxb) {
			for (int kk = 0; kk <= nT(); kk++) {
				delete [] _hidxb[kk];
			}
			delete [] _hidxb;
			}	
		}

		/**
		 * \brief Copy constructor
		 *
		 * Copying is not allowed.
		 */
		HPMPCSolver(HPMPCSolver const&) = delete;

		/**
		 * \brief Move constructor
		 *
		 * Move-construction is ok.
		 */
		HPMPCSolver(HPMPCSolver&& rhs)
		:	_nx(std::move(rhs._nx))
		,	_nu(std::move(rhs._nu))
		,	_nb(std::move(rhs._nb))
		,	_ng(std::move(rhs._ng))
		,	_workspace(std::move(rhs._workspace))
		,	_stat(std::move(rhs._stat))
		,	_mu0(rhs._mu0)
		,	_muTol(rhs._muTol)
		,	_warmStart(rhs._warmStart)
		,	infinity_(rhs.infinity_)
		,	stageBounds_(std::move(rhs.stageBounds_))
		,	terminalStageBounds_(rhs.terminalStageBounds_)
		,	lb_(std::move(rhs.lb_))
		,	ub_(std::move(rhs.ub_))
		,	lg_(std::move(rhs.lg_))
		,	ug_(std::move(rhs.ug_))
		{
			lb_.back() = terminalStageBounds_.lb.data();
			ub_.back() = terminalStageBounds_.ub.data();
			lg_.back() = terminalStageBounds_.lg.data();
			ug_.back() = terminalStageBounds_.ug.data();
			
			_hidxb = rhs._hidxb;
			rhs._hidxb = nullptr;
		}

		HPMPCSolver& operator= (HPMPCSolver const&) = delete;
		HPMPCSolver& operator= (HPMPCSolver &&) = delete;

		void Solve(Problem const& p, Solution& s)
		{
			if (p.nT() != nT())
				throw std::invalid_argument("HPMPCSolver::Solve(): problem has different number of stages than solver.");

			// Copy problem bounds and replace infinities.
			for (std::size_t i = 0; i < nT(); ++i)
			{
				CopyReplaceInfinity(p.lb_data()[i], stageBounds_[i].lb);
				CopyReplaceInfinity(p.ub_data()[i], stageBounds_[i].ub);
				CopyReplaceInfinity(p.lg_data()[i], stageBounds_[i].lg);
				CopyReplaceInfinity(p.ug_data()[i], stageBounds_[i].ug);
			}

			CopyReplaceInfinity(p.lb_data()[nT()], terminalStageBounds_.lb);
			CopyReplaceInfinity(p.ub_data()[nT()], terminalStageBounds_.ub);
			CopyReplaceInfinity(p.lg_data()[nT()], terminalStageBounds_.lg);
			CopyReplaceInfinity(p.ug_data()[nT()], terminalStageBounds_.ug);

			// Call HPMPC
			int num_iter = 0;

			auto const ret = hpmpc_wrapper::c_order_d_ip_ocp_hard_tv(&num_iter, getMaxIter(), _mu0, _muTol, nT(),
					_nx.data(), _nu.data(), _nb.data(), _hidxb, _ng.data(), nT(), _warmStart ? 1 : 0, p.A_data(), p.B_data(), p.b_data(),
					p.Q_data(), p.S_data(), p.R_data(), p.q_data(), p.r_data(), lb_.data(), ub_.data(), p.C_data(), p.D_data(),
					lg_.data(), ug_.data(), s.x_data(), s.u_data(), s.pi_data(), s.lam_data(), s.t_data(), s.inf_norm_res_data(),
					_workspace.data(), _stat[0].data());
			if (ret != 0)
			{
				throw HpmpcUnsolvedQpException(p, ret);
			}

			s.setNumIter(num_iter);

			// Warmstarting disabled on purpose.
			// On AMD K8 (hpmpc compiled for SSE3), WITHOUT warmstarting it is significantly
			// FASTER (9ms vs 14ms per time step) than with warmstarting. I am curious why.
			_warmStart = false;
		}

		std::size_t nT() const noexcept { return _nx.size() - 1; }
		std::size_t getMaxIter() const noexcept { return _stat.size(); }

		double getMuTol() const noexcept { return _muTol; }

		void setMuTol(double val)
		{
			if (val <= 0.)
				throw std::invalid_argument("mu tolerance for hpmpc must be positive");

			_muTol = val;
		}

		/**
		 * \brief Get value which is used to replace infinities in problem bounds to make HPMPC happy.
		 */
		Scalar getInfinity() const
		{
			return infinity_;
		}

		/**
		 * \brief Set value which is used to replace infinities in problem bounds to make HPMPC happy.
		 *
		 * Normal infinity will not, so don't forget to change it to something else if you are going to use infinite bounds!
		 */
		void setInfinity(Scalar val)
		{
			if (!(val > 0.))
				throw std::invalid_argument("HPMPCSolver::setInvinity(): infinity must be greater than 0.");

			infinity_ = val;
		}

	private:
		template <std::size_t NB, std::size_t NG>
		struct Bounds
		{
			std::array<Scalar, NB> lb;
			std::array<Scalar, NB> ub;
			std::array<Scalar, NG> lg;
			std::array<Scalar, NG> ug;
		};

		typedef Bounds<NX + NU, NC> StageBounds;
		typedef Bounds<NX, NCT> TerminalStageBounds;

		template <std::size_t N>
		void CopyReplaceInfinity(Scalar const * src, std::array<Scalar, N>& dst)
		{
			auto const inf = infinity_;
			std::transform(src, src + N, dst.begin(),
					[inf] (Scalar x) { return std::isinf(x) ? (x < 0. ? -inf : inf) : x; });
		}

		// Array of NX sizes
		std::vector<int> _nx;

		// Array of NU sizes
		std::vector<int> _nu;

		// Array of NB (bound constraints) sizes
		std::vector<int> _nb;

		// Array of NG (path constraints) sizes
		std::vector<int> _ng;

		// Workspace for HPMPC functions
		std::vector<char> _workspace;

		// Iteration statistics. HPMPC returns 5 double numbers per iteration.
		typedef std::array<double, 5> IterStat;
		std::vector<IterStat> _stat;

		double _mu0 = 0.;
		double _muTol = 1e-10;
		bool _warmStart = false;

		// Value which is used to replace infinities in problem bounds to make HPMPC happy
		Scalar infinity_ = std::numeric_limits<Scalar>::infinity();

		// Bounds with infinities replaced to be passed to HPMPC
		std::vector<StageBounds> stageBounds_;

		// Terminal stage bounds with infinities replaced to be passed to HPMPC
		TerminalStageBounds terminalStageBounds_;

		// Lower state bound pointers passed to HPMPC
		std::vector<Scalar const *> lb_;

		// Upper state bound pointers passed to HPMPC
		std::vector<Scalar const *> ub_;

		// Lower constraint bound pointers passed to HPMPC
		std::vector<Scalar const *> lg_;

		// Upper constraint bound pointers passed to HPMPC
		std::vector<Scalar const *> ug_;
		
		// Additional data for new hpmpc interface
		int** _hidxb; 
		
	};
}
