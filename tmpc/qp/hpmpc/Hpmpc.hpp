#pragma once

#include <tmpc/qp/QpSolverException.hpp>

#include <c_interface.h>

#include <cstdlib>
#include <cmath>


namespace tmpc :: hpmpc
{
    namespace detail
    {
        /// @brief Check if all elements in x[0], x[1], ..., x[n-1] are finite.
        template <typename Real>
        inline bool isfinite(Real const * x, size_t n)
        {
            size_t i = 0;
            while (i < n && std::isfinite(x[i]))
                ++i;
                
            return i == n;
        }
    }


    inline void c_order_ip_ocp_hard_tv(
        int& kk, int k_max, double mu0, double mu_tol,
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
        // Check arguments
        for (int i = 0; i < N; ++i)
        {
            if (!detail::isfinite(A[i], nx[i + 1] * nx[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("A is inf or nan"));

            if (!detail::isfinite(B[i], nx[i + 1] * nu[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("B is inf or nan"));

            if (!detail::isfinite(b[i], nx[i + 1]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("b is inf or nan"));
        }

        for (int i = 0; i <= N; ++i)
        {
            if (!detail::isfinite(Q[i], nx[i] * nx[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("Q is inf or nan"));

            if (!detail::isfinite(S[i], nu[i] * nx[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("S is inf or nan"));

            if (!detail::isfinite(R[i], nu[i] * nu[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("R is inf or nan"));

            if (!detail::isfinite(q[i], nx[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("q is inf or nan"));

            if (!detail::isfinite(r[i], nu[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("r is inf or nan"));

            if (!detail::isfinite(lb[i], nb[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("lb is inf or nan"));

            if (!detail::isfinite(ub[i], nb[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("ub is inf or nan"));

            if (!detail::isfinite(C[i], ng[i] * nx[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("C is inf or nan"));

            if (!detail::isfinite(D[i], ng[i] * nu[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("D is inf or nan"));

            if (!detail::isfinite(lg[i], ng[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("lg is inf or nan"));

            if (!detail::isfinite(ug[i], ng[i]))
                TMPC_THROW_EXCEPTION(std::invalid_argument("ug is inf or nan"));
        }

        auto const ret = ::c_order_d_ip_ocp_hard_tv(
            &kk, k_max, mu0, mu_tol,
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

        // Check that the number of iterations doex not exceed the maximum.
        // assert(kk <= k_max);

        // HPMPC returns garbage in the kk variable.
        // We set it to k_max as a workaround.
        kk = k_max;

        if (ret != 0)
            TMPC_THROW_EXCEPTION(QpSolverException {}
                << HpmpcReturnCodeInfo {ret}
                << boost::errinfo_api_function("c_order_d_ip_ocp_hard_tv")
            );
    }

    
    template <typename Real>
    int ip_ocp_hard_tv_work_space_size_bytes(int N, int const *nx, int const *nu, int const *nb, int const * const * hidxb, int const *ng, int N2);


    template <>
    inline int ip_ocp_hard_tv_work_space_size_bytes<double>(int N, int const *nx, int const *nu, int const *nb, int const * const * hidxb, int const *ng, int N2)
    {
        return ::hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(
            N, const_cast<int*>(nx), const_cast<int*>(nu),  const_cast<int*>(nb), const_cast<int **>(hidxb), const_cast<int*>(ng), N2);
    }
}
