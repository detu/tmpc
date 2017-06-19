/*
 * HPMPCSolver.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: kotlyar
 */

#include <tmpc/qp/HpmpcWorkspace.hpp>

#include "HPMPCProblemExport.hpp"

#include <hpmpc/c_interface.h>

#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>
#include <iomanip>

namespace tmpc
{
	HpmpcException::HpmpcException(int code)
	:	QpSolverException("HPMPC"),
		_code(code),
		msg_(std::string(QpSolverException::what()) + "\nReturn code " + std::to_string(code))
	{
	}

	static int ip_ocp_hard_tv(
		int *kk, int k_max, double mu0, double mu_tol,
		int N, int const *nx, int const *nu, int const *nb, int const *ng,
		int warm_start,
		double const * const *A, double const * const *B, double const * const *b,
		double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r,
		double const * const *lb, double const * const *ub,
		double const * const *C, double const * const *D, double const * const *lg, double const * const *ug,
		double * const *x, double * const *u, double * const *pi, double * const *lam, double * const *t,
		double *inf_norm_res,
		void *work0,
		double *stat)
	{
		return ::c_order_d_ip_ocp_hard_tv(
			kk, k_max, mu0, mu_tol,
			N, nx, nu, nb, ng,
			warm_start,
			A, B, b,
			Q, S, R, q, r,
			lb, ub,
			C, D, lg, ug,
			x, u, pi, lam, t,
			inf_norm_res,
			work0,
			stat);
	}

	/*
	static int ip_ocp_hard_tv(
		int *kk, int k_max, float mu0, float mu_tol,
		int N, int const *nx, int const *nu, int const *nb, int const *ng,
		int warm_start,
		float const * const *A, float const * const *B, float const * const *b,
		float const * const *Q, float const * const *S, float const * const *R, float const * const *q, float const * const *r,
		float const * const *lb, float const * const *ub,
		float const * const *C, float const * const *D, float const * const *lg, float const * const *ug,
		float * const *x, float * const *u, float * const *pi, float * const *lam, float * const *t,
		float *inf_norm_res,
		void *work0,
		float *stat)
	{
		::c_order_s_ip_ocp_hard_tv(
			kk, k_max, mu0, mu_tol,
			N, nx, nu, nb, ng,
			warm_start,
			A, B, b,
			Q, S, R, q, r,
			lb, ub,
			C, D, lg, ug,
			x, u, pi, lam, t,
			inf_norm_res,
			work0,
			stat);
	}
	*/

	int fortran_order_d_ip_ocp_hard_tv(int *kk, int k_max, double mu0, double mu_tol, int N, int *nx, int *nu, int *nb, int *ng, int warm_start, double **A, double **B, double **b, double **Q, double **S, double **R, double **q, double **r, double **lb, double **ub, double **C, double **D, double **lg, double **ug, double **x, double **u, double **pi, double **lam, double **t, double *inf_norm_res, void *work0, double *stat);

	template <typename Scalar>
	static int ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int const *ng);

	template <>
	int ip_ocp_hard_tv_work_space_size_bytes<double>(int N, int const *nx, int const *nu, int const *nb, int const *ng)
	{
		return ::hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(N, nx, nu, nb, ng);
	}

	template <typename Scalar_>
	void HpmpcWorkspace<Scalar_>::solve()
	{
		if (size() > 0)
		{
			// Number of QP steps for HPMPC
			auto const N = size() - 1;

			// Make sure we have enough workspace.
			solverWorkspace_.resize(ip_ocp_hard_tv_work_space_size_bytes<Scalar_>(
					static_cast<int>(N), nx_.data(), nu_.data(), nb_.data(), ng_.data()));

			// Call HPMPC
			auto const ret = ip_ocp_hard_tv(&numIter_, maxIter(), mu_, muTol_, N,
					nx_.data(), nu_.data(), nb_.data(), ng_.data(), _warmStart ? 1 : 0, A_.data(), B_.data(), b_.data(),
					Q_.data(), S_.data(), R_.data(), q_.data(), r_.data(), lb_.data(), ub_.data(), C_.data(), D_.data(),
					lg_.data(), ug_.data(), x_.data(), u_.data(), pi_.data(), lam_.data(), t_.data(), infNormRes_.data(),
					solverWorkspace_.data(), stat_[0].data());

			if (ret != 0)
			{
				throw HpmpcException(ret);
			}
		}
	}

	template <typename Scalar_>
	HpmpcWorkspace<Scalar_>::Stage::Stage(QpSize const& sz, size_t nx_next)
	:	size_(sz)
	// Initialize all numeric data to NaN so that if an uninitialized object
	// by mistake used in calculations is easier to detect.
	,	Q_(sz.nx(), sz.nx(), sNaN())
	,	R_(sz.nu(), sz.nu(), sNaN())
	,	S_(sz.nx(), sz.nu(), sNaN())
	,	q_(sz.nx(), sNaN())
	,	r_(sz.nu(), sNaN())
	,	A_(nx_next, sz.nx(), sNaN())
	,	B_(nx_next, sz.nu(), sNaN())
	,	b_(nx_next, sNaN())
	,	C_(sz.nc(), sz.nx(), sNaN())
	,	D_(sz.nc(), sz.nu(), sNaN())
	,	lb_(sz.nu() + sz.nx(), sNaN())
	,	ub_(sz.nu() + sz.nx(), sNaN())
	,	lbd_(sz.nc(), sNaN())
	,	ubd_(sz.nc(), sNaN())
	,	x_(sz.nx(), sNaN())
	,	u_(sz.nu(), sNaN())
	,	pi_(nx_next, sNaN())
	,	lam_(2 * sz.nc() + 2 * (sz.nx() + sz.nu()), sNaN())
	,	t_(2 * sz.nc() + 2 * (sz.nx() + sz.nu()), sNaN())
	{
	}

	template <typename Scalar_>
	void HpmpcWorkspace<Scalar_>::addStage(QpSize const& sz, size_t nx_next)
	{
		stage_.emplace_back(sz, nx_next);
		auto& st = stage_.back();

		nx_.push_back(sz.nx());
		nu_.push_back(sz.nu());
		nb_.push_back(sz.nx() + sz.nu());
		ng_.push_back(sz.nc());

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
		t_.push_back(st.t_data());
	}

	template <typename Scalar_>
	void HpmpcWorkspace<Scalar_>::preallocateStages(size_t nt)
	{
		stage_.reserve(nt);

		nx_.reserve(nt);
		nu_.reserve(nt);
		nb_.reserve(nt);
		ng_.reserve(nt);

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
		t_.reserve(nt);
	}

	// Explicit instantiation of HpmpcWorkspace for supported data types and storage orders.
	template class HpmpcWorkspace<double>;
}
