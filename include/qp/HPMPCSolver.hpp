#pragma once

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

	template<unsigned NX_, unsigned NU_, unsigned NC_, unsigned NCT_>
	class HPMPCSolver
	{
	public:
		static unsigned const NX = NX_;
		static unsigned const NU = NU_;
		static unsigned const NZ = NX + NU;
		static unsigned const NC = NC_;
		static unsigned const NCT = NCT_;

		typedef HPMPCProblem<NX, NU, NC, NCT> Problem;
		typedef HPMPCSolution<NX, NU, NC, NCT> Solution;

		HPMPCSolver(std::size_t nt, int max_iter = 100)
		:	_nx(nt + 1, NX)
		,	_nu(nt + 1, NU)
		,	_nb(nt + 1, NU + NX)
		,	_ng(nt + 1, NC)
		,	_stat(max_iter)
		{
			//_nx.front() = 0;	// HPMPC wants nx[0] = 0 for MPC problems (seems that it does not support lbx[0] == ubx[0])
			_nu.back() = 0;
			_nb.back() = NX;
			_ng.back() = NCT;

			// Allocate workspace
			_workspace.resize(hpmpc_wrapper::d_ip_ocp_hard_tv_work_space_size_bytes(
					static_cast<int>(nt), _nx.data(), _nu.data(), _nb.data(), _ng.data()));
		}

		void Solve(Problem const& p, Solution& s)
		{
			int num_iter = 0;

			bool const x0_equality = p.is_x0_equality_constrained();
			_nx[0] = x0_equality ? 0 : NX;
			_nb[0] = _nx[0] + _nu[0];

			hpmpc_wrapper::c_order_d_ip_ocp_hard_tv(&num_iter, getMaxIter(), _mu0, _muTol, nT(),
					_nx.data(), _nu.data(), _nb.data(), _ng.data(), _warmStart ? 1 : 0, p.A_data(), p.B_data(), p.b_data(),
					p.Q_data(), p.S_data(), p.R_data(), p.q_data(), p.r_data(), p.lb_data(), p.ub_data(), p.C_data(), p.D_data(),
					p.lg_data(), p.ug_data(), s.x_data(), s.u_data(), s.pi_data(), s.lam_data(), s.t_data(), s.inf_norm_res_data(),
					_workspace.data(), _stat[0].data());

			if (x0_equality)
				s.set_x(0, p.get_x_min(0));

			// Warmstarting disabled on purpose.
			// 1. After the Simuling model is executed about 3 times, next runs produce
			// the "HPMPC returned -1" error. This happens randomly. To make sure that
			// this bug has nothing to do with warmstarting, I disable it.
			//_warmStart = true;
			// 2.On AMD K8 (hpmpc compiled for SSE3), WITHOUT warmstarting it is significantly
			// FASTER (9ms vs 14ms per time step) than with warmstarting. I am curious why.
			_warmStart = false;
		}

		std::size_t nT() const noexcept { return _nx.size() - 1; }
		std::size_t getMaxIter() const noexcept { return _stat.size(); }

	private:
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
	};
}
