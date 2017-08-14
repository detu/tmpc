#include <tmpc/qp/HpipmWorkspace.hpp>

#include <stdexcept>
#include <sstream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <vector>

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

	template <typename Real_>
	void HpipmWorkspace<Real_>::solve()
	{
		if (stage_.size() > 1)
		{
			// Number of QP steps for HPIPM
			auto const N = stage_.size() - 1;

			// Convert the problem
			HPIPM::cvt_colmaj_to_ocp_qp(
				A_.data(), B_.data(), b_.data(), 
				Q_.data(), S_.data(), R_.data(), q_.data(), r_.data(), 
				hidxb_.data(), lb_.data(), ub_.data(), 
				C_.data(), D_.data(), lg_.data(), ug_.data(), &ocpQp_);

			// Call HPIPM
			auto const ret = 0;
			HPIPM::solve_ipm_hard_ocp_qp(&ocpQp_, &ocpQpSol_, &solverWorkspace_);

			if (ret != 0)
			{
				throw HpipmException(ret);
			}

			// Convert the solution
			HPIPM::cvt_ocp_qp_sol_to_colmaj(&ocpQp_, &ocpQpSol_, 
				u_.data(), x_.data(), pi_.data(), lam_lb_.data(), lam_ub_.data(), lam_lg_.data(), lam_ug_.data());
		}
	}

	template <typename Real_>
	HpipmWorkspace<Real_>::Stage::Stage(QpSize const& sz, size_t nx_next)
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
	void HpipmWorkspace<Real_>::addStage(QpSize const& sz, size_t nx_next)
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
		lam_lb_.push_back(st.lam_lb_data());
		lam_ub_.push_back(st.lam_ub_data());
		lam_lg_.push_back(st.lam_lg_data());
		lam_ug_.push_back(st.lam_ug_data());
	}

	template <typename Real_>
	void HpipmWorkspace<Real_>::preallocateStages(size_t nt)
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
		lam_lb_.reserve(nt);
		lam_ub_.reserve(nt);
		lam_lg_.reserve(nt);
		lam_ug_.reserve(nt);
	}

	// Explicit instantiation of HpipmWorkspace for supported data types and storage orders.
	template class HpipmWorkspace<double>;
}
