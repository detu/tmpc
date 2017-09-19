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

	template <>
	int Hpmpc<double>::c_order_ip_ocp_hard_tv(
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

	template <>
	int Hpmpc<double>::ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int const * const * hidxb, int const *ng, int N2)
	{
		return ::hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(
			N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), const_cast<int **>(hidxb), const_cast<int*>(ng), N2);
	}
}
