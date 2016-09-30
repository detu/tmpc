#pragma once

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>

#include "HPMPCProblem.hpp"
#include "HPMPCSolution.hpp"

namespace tmpc
{
	namespace hpmpc_wrapper
	{
		void c_order_d_ip_ocp_hard_tv(	int *kk, int k_max, double mu0, double mu_tol,
										int N, int const *nx, int const *nu, int const *nb, int const *ng,
										int warm_start,
										double const * const *A, double const * const *B, double const * const *b,
										double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r,
										double const * const *lb, double const * const *ub,
										double const * const *C, double const * const *D, double const * const *lg, double const * const *ug,
										double * const *x, double * const *u, double * const *pi, double * const *lam, double * const *t,
										double *inf_norm_res,
										void *work0,
										double *stat );

		int d_ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int const *ng);
	}

	/**
	 * \brief Multistage QP solver using qpOASES
	 *
	 * \tparam K a class implementing the Kernel concept
	 */
	template <typename K>
	class HPMPCSolver
	{
	public:
		static auto constexpr NX = K::NX;
		static auto constexpr NU = K::NU;
		static auto constexpr NZ = K::NX + K::NU;
		static auto constexpr NC = K::NC;
		static auto constexpr NCT = K::NCT;

		typedef typename K::Scalar Scalar;

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

			// Allocate workspace
			_workspace.resize(hpmpc_wrapper::d_ip_ocp_hard_tv_work_space_size_bytes(
					static_cast<int>(nt), _nx.data(), _nu.data(), _nb.data(), _ng.data()));
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
		HPMPCSolver(HPMPCSolver&& rhs) = default;

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

			hpmpc_wrapper::c_order_d_ip_ocp_hard_tv(&num_iter, getMaxIter(), _mu0, _muTol, nT(),
					_nx.data(), _nu.data(), _nb.data(), _ng.data(), _warmStart ? 1 : 0, p.A_data(), p.B_data(), p.b_data(),
					p.Q_data(), p.S_data(), p.R_data(), p.q_data(), p.r_data(), lb_.data(), ub_.data(), p.C_data(), p.D_data(),
					lg_.data(), ug_.data(), s.x_data(), s.u_data(), s.pi_data(), s.lam_data(), s.t_data(), s.inf_norm_res_data(),
					_workspace.data(), _stat[0].data());

			// Warmstarting disabled on purpose.
			// 1. After the Simulink model is executed about 3 times, next runs produce
			// the "HPMPC returned -1" error. This happens randomly. To make sure that
			// this bug has nothing to do with warmstarting, I disable it.
			//_warmStart = true;
			// 2. On AMD K8 (hpmpc compiled for SSE3), WITHOUT warmstarting it is significantly
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
	};
}
