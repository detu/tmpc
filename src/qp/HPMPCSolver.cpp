/*
 * HPMPCSolver.cpp
 *
 *  Created on: Jun 20, 2016
 *      Author: kotlyar
 */

#include "HPMPCProblemExport.hpp"

#include <hpmpc/c_interface.h>

#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>
#include <iomanip>

namespace tmpc
{
	namespace hpmpc_wrapper
	{
		void throw_hpmpc_error(int err_code, std::string const& s = "")
		{
			std::ostringstream msg;
			msg << "HPMPC error: return code = " << err_code << " . " << s << ".";
			throw std::runtime_error(msg.str());
		}

		int c_order_d_ip_ocp_hard_tv(
									int *kk, int k_max, double mu0, double mu_tol,
									int N, int const *nx, int const *nu, int const *nb, int **hidxb, int const *ng, int N2,
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
			//int rr = feenableexcept(FE_ALL_EXCEPT);
			//double nan = std::numeric_limits<double>::signaling_NaN();
			//double xx = nan + 1;
			int const ret = ::c_order_d_ip_ocp_hard_tv(
					kk, k_max, mu0, mu_tol,
					N, const_cast<int*>(nx), const_cast<int*>(nu), const_cast<int*>(nb), hidxb, const_cast<int*>(ng), N2,
					warm_start,
					const_cast<double**>(A), const_cast<double**>(B), const_cast<double**>(b),
					const_cast<double**>(Q), const_cast<double**>(S), const_cast<double**>(R), const_cast<double**>(q), const_cast<double**>(r),
					const_cast<double**>(lb), const_cast<double**>(ub),
					const_cast<double**>(C), const_cast<double**>(D), const_cast<double**>(lg), const_cast<double**>(ug),
					const_cast<double**>(x), const_cast<double**>(u), const_cast<double**>(pi), const_cast<double**>(lam), 
					inf_norm_res,
					work0,
					stat);
			
			//int const ret = 1;
			if (ret != 0)
			{
				using namespace hpmpc_problem_export;

				{
					std::ofstream os("failed_qp_hpmpc.m");
					os << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);

					MATLABFormatter f(os, "qp.");

					print_c_order_d_ip_ocp_hard_tv(f,
						k_max, mu0, mu_tol,
						N, nx, nu, nb, ng,
						warm_start,
						A, B, b,
						Q, S, R, q, r,
						lb, ub,
						C, D, lg, ug, x, u);
				}

				{
					std::ofstream os("failed_qp_hpmpc.c");
					os << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);

					CFormatter f(os);

					print_c_order_d_ip_ocp_hard_tv(f,
						k_max, mu0, mu_tol,
						N, nx, nu, nb, ng,
						warm_start,
						A, B, b,
						Q, S, R, q, r,
						lb, ub,
						C, D, lg, ug, x, u);
				}

				//throw_hpmpc_error(ret, "The QP dumped to failed_qp_hpmpc.m and failed_qp_hpmpc.c");
			}

			return ret;
		}

		int fortran_order_d_ip_ocp_hard_tv(int *kk, int k_max, double mu0, double mu_tol, int N, int *nx, int *nu, int *nb, int *ng, int warm_start, double **A, double **B, double **b, double **Q, double **S, double **R, double **q, double **r, double **lb, double **ub, double **C, double **D, double **lg, double **ug, double **x, double **u, double **pi, double **lam, double **t, double *inf_norm_res, void *work0, double *stat);

		int d_ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int **hidxb, int const *ng, int N2)
		{
			return ::hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), hidxb, const_cast<int*>(ng), N2);
		}
	}
}
