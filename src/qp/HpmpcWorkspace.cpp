/*
 * HPMPCSolver.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: kotlyar
 */

#include <tmpc/qp/HpmpcWorkspace.hpp>

#include "HPMPCProblemExport.hpp"

#include <c_interface.h>

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
		int N, int const *nx, int const *nu, int const *nb, int const * const *hidxb, int const *ng, int N2,
		int warm_start,
		double const * const *A, double const * const *B, double const * const *b,
		double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r,
		double const * const *lb, double const * const *ub,
		double const * const *C, double const * const *D, double const * const *lg, double const * const *ug,
		double * const *x, double * const *u, double * const *pi, double * const *lam,
		double *inf_norm_res,
		void *work0,
		double *stat)
	{
		return ::c_order_d_ip_ocp_hard_tv(
			kk, k_max, mu0, mu_tol,
			N, const_cast<int*>(nx), const_cast<int*>(nu), const_cast<int*>(nb), const_cast<int **>(hidxb), const_cast<int*>(ng), N2,
			warm_start,
			const_cast<double**>(A), const_cast<double**>(B), const_cast<double**>(b),
			const_cast<double**>(Q), const_cast<double**>(S), const_cast<double**>(R), const_cast<double**>(q), const_cast<double**>(r),
			const_cast<double**>(lb), const_cast<double**>(ub),
			const_cast<double**>(C), const_cast<double**>(D), const_cast<double**>(lg), const_cast<double**>(ug),
			const_cast<double**>(x), const_cast<double**>(u), const_cast<double**>(pi), const_cast<double**>(lam), 
			inf_norm_res,
			work0,
			stat);
	}

	template <typename Real>
	static int ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int const * const * hidxb, int const *ng, int N2);

	template <>
	int ip_ocp_hard_tv_work_space_size_bytes<double>(int N, int const *nx, int const *nu, int const *nb, int const * const * hidxb, int const *ng, int N2)
	{
		return ::hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(
			N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), const_cast<int **>(hidxb), const_cast<int*>(ng), N2);
	}

	template <typename Real_>
	void HpmpcWorkspace<Real_>::solve()
	{
		if (stage_.size() > 0)
		{
			// Number of QP steps for HPMPC
			auto const N = stage_.size() - 1;

			// Make sure we have enough workspace.
			solverWorkspace_.resize(ip_ocp_hard_tv_work_space_size_bytes<Real_>(
					static_cast<int>(N), nx_.data(), nu_.data(), nb_.data(), hidxb_.data(), ng_.data(), static_cast<int>(N)));

			// Call HPMPC
			auto const ret = ip_ocp_hard_tv(&numIter_, maxIter(), mu_, muTol_, N,
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

	template <typename Real_>
	HpmpcWorkspace<Real_>::Stage::Stage(QpSize const& sz, size_t nx_next)
	:	size_(sz)
	,	hidxb_(sz.nu() + sz.nx())
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
	{
		int n = 0;
		std::generate(hidxb_.begin(), hidxb_.end(), [&n] { return n++; });
	}

	template <typename Real_>
	void HpmpcWorkspace<Real_>::addStage(QpSize const& sz, size_t nx_next)
	{
		stage_.emplace_back(sz, nx_next);
		auto& st = stage_.back();

		nx_.push_back(sz.nx());
		nu_.push_back(sz.nu());
		nb_.push_back(sz.nx() + sz.nu());
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

	template <typename Real_>
	void HpmpcWorkspace<Real_>::preallocateStages(size_t nt)
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

	// Explicit instantiation of HpmpcWorkspace for supported data types and storage orders.
	template class HpmpcWorkspace<double>;
}
