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
	int Hpipm<double>::memsize_ocp_qp(int N, int const *nx, int const *nu, int const *nb, int const *ng, int const * ns)
	{
		return ::d_memsize_ocp_qp(
			N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), const_cast<int*>(ng), const_cast<int*>(ns));
	}

	void Hpipm<double>::create_ocp_qp(int N, int const *nx, int const *nu, int const *nb, int const *ng, int const * ns, ocp_qp *qp, void *memory)
	{
		::d_create_ocp_qp(N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), const_cast<int*>(ng), const_cast<int*>(ns), qp, memory);
	}

	void Hpipm<double>::cvt_colmaj_to_ocp_qp(
		double const * const *A, double const * const *B, double const * const *b, 
		double const * const *Q, double const * const *S, double const * const *R, double const * const *q, double const * const *r, 
		int const * const *idxb, double const * const *lb, double const * const *ub, 
		double const * const *C, double const * const *D, double const * const *lg, double const * const *ug, 
		double const * const *Zl, double const * const *Zu, double const * const *zl, double const * const *zu, int const * const *idxs,
		ocp_qp * qp)
	{
		::d_cvt_colmaj_to_ocp_qp(
			const_cast<double **>(A), const_cast<double **>(B), const_cast<double **>(b), 
			const_cast<double **>(Q), const_cast<double **>(S), const_cast<double **>(R), const_cast<double **>(q), const_cast<double **>(r), 
			const_cast<int **>(idxb), const_cast<double **>(lb), const_cast<double **>(ub), 
			const_cast<double **>(C), const_cast<double **>(D), const_cast<double **>(lg), const_cast<double **>(ug), 
			const_cast<double **>(Zl), const_cast<double **>(Zu), const_cast<double **>(zl), const_cast<double **>(zu), const_cast<int **>(idxs), 
			qp
		);
	}

	int Hpipm<double>::memsize_ocp_qp_sol(int N, int const *nx, int const *nu, int const *nb, int const *ng, int const *ns)
	{
		return ::d_memsize_ocp_qp_sol(N, const_cast<int *>(nx), const_cast<int *>(nu), const_cast<int *>(nb), const_cast<int *>(ng), const_cast<int *>(ns));
	}

	void Hpipm<double>::create_ocp_qp_sol(int N, int const *nx, int const *nu, int const *nb, int const *ng, int const *ns, ocp_qp_sol *qp_sol, void *memory)
	{
		::d_create_ocp_qp_sol(N, const_cast<int *>(nx), const_cast<int *>(nu), const_cast<int *>(nb), const_cast<int *>(ng), const_cast<int *>(ns), qp_sol, memory);
	}

	void Hpipm<double>::cvt_ocp_qp_sol_to_colmaj(ocp_qp const *qp, ocp_qp_sol const *qp_sol, 
		double * const *u, double * const *x, double * const * ls, double * const * us,
		double * const *pi, 
		double * const *lam_lb, double * const *lam_ub, double * const *lam_lg, double * const *lam_ug,
		double * const *lam_ls, double * const *lam_us)
	{
		::d_cvt_ocp_qp_sol_to_colmaj(const_cast<ocp_qp *>(qp), const_cast<ocp_qp_sol *>(qp_sol), 
			const_cast<double **>(u), const_cast<double **>(x), const_cast<double **>(ls), const_cast<double **>(us), 
			const_cast<double **>(pi), 
			const_cast<double **>(lam_lb), const_cast<double **>(lam_ub), const_cast<double **>(lam_lg), const_cast<double **>(lam_ug),
			const_cast<double **>(lam_ls), const_cast<double **>(lam_us)
		);
	}


	int Hpipm<double>::memsize_ocp_qp_ipm(ocp_qp const *qp, ocp_qp_ipm_arg const *arg)
	{
		return ::d_memsize_ocp_qp_ipm(const_cast<ocp_qp *>(qp), const_cast<ocp_qp_ipm_arg *>(arg));
	}


	void Hpipm<double>::create_ocp_qp_ipm(ocp_qp const *qp, ocp_qp_ipm_arg const *arg, ocp_qp_ipm_workspace *ws, void *mem)
	{
		::d_create_ocp_qp_ipm(const_cast<ocp_qp *>(qp), const_cast<ocp_qp_ipm_arg *>(arg), ws, mem);
	}


	int Hpipm<double>::solve_ocp_qp_ipm(ocp_qp const *qp, ocp_qp_sol *qp_sol, ocp_qp_ipm_arg const *arg, ocp_qp_ipm_workspace *ws)
	{
		return ::d_solve_ocp_qp_ipm(const_cast<ocp_qp *>(qp), qp_sol, const_cast<ocp_qp_ipm_arg *>(arg), ws);
	}


	HpipmException::HpipmException(int code)
	:	QpSolverException("HPIPM"),
		_code(code),
		msg_(std::string(QpSolverException::what()) + "\nReturn code " + std::to_string(code))
	{
	}
}
