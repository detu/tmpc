#include <tmpc/qp/HpipmWorkspace.hpp>

#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <vector>
#include <algorithm>

namespace tmpc
{
	int Hpipm<double>::memsize_ocp_qp(int N, int const *nx, int const *nu, int const *nb, int const *ng)
	{
		return ::d_memsize_ocp_qp(
			N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), const_cast<int*>(ng));
	}

	void Hpipm<double>::create_ocp_qp(int N, int const *nx, int const *nu, int const *nb, int const *ng, ocp_qp *qp, void *memory)
	{
		::d_create_ocp_qp(N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), const_cast<int*>(ng), qp, memory);
	}

	void Hpipm<double>::cvt_colmaj_to_ocp_qp(
		double const * const *A, double const * const *B, double const * const *b, 
		double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r, 
		int const * const *idxb, double const * const *lb, double const * const *ub, 
		double const * const *C, double const * const *D, double const * const *lg, double const * const *ug, ocp_qp * qp)
	{
		::d_cvt_colmaj_to_ocp_qp(
			const_cast<double **>(A), const_cast<double **>(B), const_cast<double **>(b), 
			const_cast<double **>(Q), const_cast<double **>(S), const_cast<double **>(R), const_cast<double **>(q), const_cast<double **>(r), 
			const_cast<int **>(idxb), const_cast<double **>(lb), const_cast<double **>(ub), 
			const_cast<double **>(C), const_cast<double **>(D), const_cast<double **>(lg), const_cast<double **>(ug), qp
		);
	}

	int Hpipm<double>::memsize_ocp_qp_sol(int N, int const *nx, int const *nu, int const *nb, int const *ng)
	{
		return ::d_memsize_ocp_qp_sol(N, const_cast<int *>(nx), const_cast<int *>(nu), const_cast<int *>(nb), const_cast<int *>(ng));
	}

	void Hpipm<double>::create_ocp_qp_sol(int N, int const *nx, int const *nu, int const *nb, int const *ng, ocp_qp_sol *qp_sol, void *memory)
	{
		::d_create_ocp_qp_sol(N, const_cast<int *>(nx), const_cast<int *>(nu), const_cast<int *>(nb), const_cast<int *>(ng), qp_sol, memory);
	}

	void Hpipm<double>::cvt_ocp_qp_sol_to_colmaj(ocp_qp const *qp, ocp_qp_sol const *qp_sol, 
		double * const *u, double * const *x, double * const *pi, double * const *lam_lb, double * const *lam_ub, double * const *lam_lg, double * const *lam_ug)
	{
		::d_cvt_ocp_qp_sol_to_colmaj(const_cast<ocp_qp *>(qp), const_cast<ocp_qp_sol *>(qp_sol), 
			const_cast<double **>(u), const_cast<double **>(x), const_cast<double **>(pi), 
			const_cast<double **>(lam_lb), const_cast<double **>(lam_ub), const_cast<double **>(lam_lg), const_cast<double **>(lam_ug));
	}

	int Hpipm<double>::memsize_ipm_hard_ocp_qp(ocp_qp const *qp, ipm_hard_ocp_qp_arg const *arg)
	{
		return ::d_memsize_ipm_hard_ocp_qp(const_cast<ocp_qp *>(qp), const_cast<ipm_hard_ocp_qp_arg *>(arg));
	}

	void Hpipm<double>::create_ipm_hard_ocp_qp(ocp_qp const *qp, ipm_hard_ocp_qp_arg const *arg, ipm_hard_ocp_qp_workspace *ws, void *mem)
	{
		::d_create_ipm_hard_ocp_qp(const_cast<ocp_qp *>(qp), const_cast<ipm_hard_ocp_qp_arg *>(arg), ws, mem);
	}

	void Hpipm<double>::solve_ipm_hard_ocp_qp(ocp_qp const *qp, ocp_qp_sol *qp_sol, ipm_hard_ocp_qp_workspace *ws)
	{
		::d_solve_ipm_hard_ocp_qp(const_cast<ocp_qp *>(qp), qp_sol, ws);
	}

	template <typename Real>
	HpipmWorkspace<Real>::OcpQp::OcpQp(int N, int const *nx, int const *nu, int const *nb, int const *ng)
	:	memory_(HPIPM::memsize_ocp_qp(N, nx, nu, nb, ng))
	{
		HPIPM::create_ocp_qp(N, nx, nu, nb, ng, &qp_, memory_.data());
	}

	template <typename Real>
	HpipmWorkspace<Real>::OcpQpSol::OcpQpSol(int N, int const *nx, int const *nu, int const *nb, int const *ng)
	:	memory_(HPIPM::memsize_ocp_qp_sol(N, nx, nu, nb, ng))
	{
		HPIPM::create_ocp_qp_sol(N, nx, nu, nb, ng, &sol_, memory_.data());
	}

	template <typename Real>
	HpipmWorkspace<Real>::IpmHardOcpQpWorkspace::IpmHardOcpQpWorkspace(typename HPIPM::ocp_qp const& qp, typename HPIPM::ipm_hard_ocp_qp_arg const& arg)
	:	memory_(HPIPM::memsize_ipm_hard_ocp_qp(&qp, &arg))
	{
		HPIPM::create_ipm_hard_ocp_qp(&qp, &arg, &workspace_, memory_.data());
	}

	HpipmException::HpipmException(int code)
	:	QpSolverException("HPIPM"),
		_code(code),
		msg_(std::string(QpSolverException::what()) + "\nReturn code " + std::to_string(code))
	{
	}
}
