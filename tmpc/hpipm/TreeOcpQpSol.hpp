#pragma once

#include <hpipm_d_tree_ocp_qp_sol.h>

#include <blasfeo/BlasfeoApi.hpp>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct TreeOcpQpSolType;


    template <>
    struct TreeOcpQpSolType<double>
    {
        using Type = ::d_tree_ocp_qp_sol;
    };
    

    template <typename Real>
    using tree_ocp_qp_sol = TreeOcpQpSolType<Real>::Type;


    inline int memsize_tree_ocp_qp_sol(d_tree_ocp_qp_dim const& dim)
    {
        return ::d_memsize_tree_ocp_qp_sol(const_cast<d_tree_ocp_qp_dim *>(&dim));
    }


    inline void create_tree_ocp_qp_sol(d_tree_ocp_qp_dim const& dim, d_tree_ocp_qp_sol& qp_sol, void * memory)
    {
        ::d_create_tree_ocp_qp_sol(const_cast<d_tree_ocp_qp_dim *>(&dim), &qp_sol, memory);
    }


    inline void cvt_tree_ocp_qp_sol_to_colmaj(d_tree_ocp_qp const& qp, d_tree_ocp_qp_sol const& qp_sol,
        double * const * u, double * const * x, double * const *ls, double * const * us,
        double * const * pi, double * const * lam_lb, double * const * lam_ub,
        double * const * lam_lg, double * const * lam_ug, double * const * lam_ls, double * const * lam_us)
    {
        ::d_cvt_tree_ocp_qp_sol_to_colmaj(const_cast<d_tree_ocp_qp *>(&qp), const_cast<d_tree_ocp_qp_sol *>(&qp_sol),
            const_cast<double **>(u), const_cast<double **>(x), const_cast<double **>(ls), const_cast<double **>(us),
            const_cast<double **>(pi), const_cast<double **>(lam_lb), const_cast<double **>(lam_ub),
            const_cast<double **>(lam_lg), const_cast<double **>(lam_ug), 
            const_cast<double **>(lam_ls), const_cast<double **>(lam_us));
    }


    inline void cvt_tree_ocp_qp_sol_to_rowmaj(d_tree_ocp_qp const& qp, d_tree_ocp_qp_sol const& qp_sol,
        double * const * u, double * const * x, double * const *ls, double * const * us,
        double * const * pi, double * const * lam_lb, double * const * lam_ub,
        double * const * lam_lg, double * const * lam_ug, double * const * lam_ls, double * const * lam_us)
    {
        ::d_cvt_tree_ocp_qp_sol_to_rowmaj(const_cast<d_tree_ocp_qp *>(&qp), const_cast<d_tree_ocp_qp_sol *>(&qp_sol),
            const_cast<double **>(u), const_cast<double **>(x), const_cast<double **>(ls), const_cast<double **>(us),
            const_cast<double **>(pi), const_cast<double **>(lam_lb), const_cast<double **>(lam_ub),
            const_cast<double **>(lam_lg), const_cast<double **>(lam_ug), 
            const_cast<double **>(lam_ls), const_cast<double **>(lam_us));
    }


    /// @brief Wrapper for tree_ocp_qp_sol struct.
    ///
    /// The implementation of the get_() functions is based on CVT_TREE_OCP_QP_SOL_TO_COLMAJ()
    /// from https://github.com/giaf/hpipm/blob/master/tree_ocp_qp/x_tree_ocp_qp_sol.c
    ///
    template <typename Real>
    struct TreeOcpQpSol
    :   tree_ocp_qp_sol<Real>
    {
        TreeOcpQpSol(d_tree_ocp_qp_dim const& dim)
        :   memory_ {new char[memsize_tree_ocp_qp_sol(dim)]}
        {
            create_tree_ocp_qp_sol(dim, *this, memory_.get());
        }


        void cvt_to_colmaj(
            d_tree_ocp_qp const& qp, d_tree_ocp_qp_sol const& qp_sol,
            Real * const * u, Real * const * x, Real * const *ls, Real * const * us,
            Real * const * pi, Real * const * lam_lb, Real * const * lam_ub,
            Real * const * lam_lg, Real * const * lam_ug, Real * const * lam_ls, Real * const * lam_us) const
        {
            cvt_tree_ocp_qp_sol_to_colmaj(qp, *this,
                u, x, ls, us, pi, lam_lb, lam_ub, lam_lg, lam_ug, lam_ls, lam_us);
        }


        void cvt_to_rowmaj(
            d_tree_ocp_qp const& qp, d_tree_ocp_qp_sol const& qp_sol,
            Real * const * u, Real * const * x, Real * const *ls, Real * const * us,
            Real * const * pi, Real * const * lam_lb, Real * const * lam_ub,
            Real * const * lam_lg, Real * const * lam_ug, Real * const * lam_ls, Real * const * lam_us) const
        {
            cvt_tree_ocp_qp_sol_to_rowmaj(qp, *this,
                u, x, ls, us, pi, lam_lb, lam_ub, lam_lg, lam_ug, lam_ls, lam_us);
        }


        void get_u(int ii, Real * vec) const
        {
            // CVT_STRVEC2VEC(nu[ii], qp_sol->ux+ii, 0, u[ii]);
            blasfeo::unpack_vec(this->dim->nu[ii], this->ux[ii], 0, vec);
        }

        
        void get_x(int ii, Real *vec) const
        {
            // CVT_STRVEC2VEC(nx[ii], qp_sol->ux+ii, nu[ii], x[ii]);
            blasfeo::unpack_vec(this->dim->nx[ii], this->ux[ii], this->dim->nu[ii], vec);
        }

        
        void get_sl(int ii, Real *vec) const
        {
            // CVT_STRVEC2VEC(ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii], ls[ii]);
            blasfeo::unpack_vec(this->dim->ns[ii], this->ux[ii], this->dim->nu[ii] + this->dim->nx[ii], vec);
        }

        
        void get_su(int ii, Real *vec) const
        {
            // CVT_STRVEC2VEC(ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii]+ns[ii], us[ii]);
            blasfeo::unpack_vec(this->dim->ns[ii], this->ux[ii], this->dim->nu[ii] + this->dim->nx[ii] + this->dim->ns[ii], vec);
        }

        
        void get_pi(int ii, Real *vec) const
        {
            // CVT_STRVEC2VEC(nx[ii+1], qp_sol->pi+ii, 0, pi[ii]);
            blasfeo::unpack_vec(this->dim->nx[ii + 1], this->pi[ii], 0, vec);
        }

        
        void get_lam_lb(int ii, Real *vec) const
        {
            // CVT_STRVEC2VEC(nb[ii], qp_sol->lam+ii, 0, lam_lb[ii]);
            blasfeo::unpack_vec(this->dim->nb[ii], this->lam[ii], 0, vec);
        }

        
        void get_lam_ub(int ii, Real *vec) const
        {
            // CVT_STRVEC2VEC(nb[ii], qp_sol->lam+ii, nb[ii]+ng[ii], lam_ub[ii]);
            blasfeo::unpack_vec(this->dim->nb[ii], this->lam[ii], this->dim->nb[ii] + this->dim->ng[ii], vec);
        }

        
        void get_lam_lg(int ii, Real *vec) const
        {
            // CVT_STRVEC2VEC(ng[ii], qp_sol->lam+ii, nb[ii], lam_lg[ii]);
            blasfeo::unpack_vec(this->dim->ng[ii], this->lam[ii], this->dim->nb[ii], vec);
        }

        
        void get_lam_ug(int ii, Real *vec) const
        {
            // CVT_STRVEC2VEC(ng[ii], qp_sol->lam+ii, 2*nb[ii]+ng[ii], lam_ug[ii]);
            blasfeo::unpack_vec(this->dim->ng[ii], this->lam[ii], 2 * this->dim->nb[ii] + this->dim->ng[ii], vec);
        }

    
    private:
        std::unique_ptr<char []> memory_;
    };
}