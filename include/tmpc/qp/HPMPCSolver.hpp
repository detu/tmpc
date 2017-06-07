#pragma once

#include "HPMPCProblem.hpp"
#include "HPMPCSolution.hpp"
#include "QpSolverException.hpp"

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace tmpc
{
	namespace hpmpc_wrapper
	{
		int c_order_d_ip_ocp_hard_tv(	int *kk, int k_max, double mu0, double mu_tol,
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

	class HpmpcUnsolvedQpException : public QpSolverException
	{
	public:
		template <typename QP>
		HpmpcUnsolvedQpException(QP const& qp, int code)
		:	QpSolverException("HPMPC", qp),
			_code(code),
			msg_(std::string(QpSolverException::what()) + "\nReturn code " + std::to_string(code))
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
	 * \tparam <Scalar_> Scalar type
	 */
	template <typename Scalar_>
	class HPMPCSolver
	{
	public:
		typedef Scalar_ Scalar;
		typedef HPMPCProblem<Scalar_> Problem;
		typedef HPMPCSolution<Scalar_> Solution;

		HPMPCSolver(int max_iter = 100)
		:	_stat(max_iter)
		{
		}

		/**
		 * \brief Takes QP problem size to preallocate workspace.
		 */
		template <typename InputIterator>
		HPMPCSolver(InputIterator sz_first, InputIterator sz_last, int max_iter = 100)
		:	_stat(max_iter)
		{
			// TODO: calculate workspace size and call
			// workspace_.reserve();
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

		HPMPCSolver& operator= (HPMPCSolver const&) = delete;
		HPMPCSolver& operator= (HPMPCSolver &&) = delete;

		void Solve(Problem const& p, Solution& s)
		{
			if (p.size() > 1)
			{
				// Number of QP steps for HPMPC
				auto const N = p.size() - 1;

				// Make sure we have enough workspace.
				_workspace.resize(hpmpc_wrapper::d_ip_ocp_hard_tv_work_space_size_bytes(
						static_cast<int>(N), p.nx_data(), p.nu_data(), p.nb_data(), p.ng_data()));

				// Call HPMPC
				int num_iter = 0;

				auto const ret = hpmpc_wrapper::c_order_d_ip_ocp_hard_tv(&num_iter, getMaxIter(), _mu0, _muTol, N,
						p.nx_data(), p.nu_data(), p.nb_data(), p.ng_data(), _warmStart ? 1 : 0, p.A_data(), p.B_data(), p.b_data(),
						p.Q_data(), p.S_data(), p.R_data(), p.q_data(), p.r_data(), p.lb_data(), p.ub_data(), p.C_data(), p.D_data(),
						p.lg_data(), p.ug_data(), s.x_data(), s.u_data(), s.pi_data(), s.lam_data(), s.t_data(), s.inf_norm_res_data(),
						_workspace.data(), _stat[0].data());

				if (ret != 0)
				{
					throw HpmpcUnsolvedQpException(p, ret);
				}

				s.setNumIter(num_iter);
			}
		}

		std::size_t getMaxIter() const noexcept { return _stat.size(); }

		double getMuTol() const noexcept { return _muTol; }

		void setMuTol(double val)
		{
			if (val <= 0.)
				throw std::invalid_argument("mu tolerance for hpmpc must be positive");

			_muTol = val;
		}
		
	private:
		// Workspace for HPMPC functions
		std::vector<char> _workspace;

		// Iteration statistics. HPMPC returns 5 double numbers per iteration.
		typedef std::array<double, 5> IterStat;
		std::vector<IterStat> _stat;

		double _mu0 = 0.;
		double _muTol = 1e-10;

		// Warmstarting disabled on purpose.
		// On AMD K8 (hpmpc compiled for SSE3), WITHOUT warmstarting it is significantly
		// FASTER (9ms vs 14ms per time step) than with warmstarting. I am curious why.
		bool _warmStart = false;
	};
}
