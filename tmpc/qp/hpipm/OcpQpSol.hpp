// **************************************************************************
//
// ocp_qp_sol API
//
// **************************************************************************
#pragma once

#include <hpipm_d_ocp_qp_sol.h>

#include <memory>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct OcpQpSolType;


    template <>
    struct OcpQpSolType<double>
    {
        using Type = ::d_ocp_qp_sol;
    };


    template <typename Real>
    using ocp_qp_sol = typename OcpQpSolType<Real>::Type;


    inline int ocp_qp_sol_memsize(d_ocp_qp_dim const& dim)
    {
        return ::d_ocp_qp_sol_memsize(const_cast<d_ocp_qp_dim *>(&dim));
    }


    inline void ocp_qp_sol_create(d_ocp_qp_dim const& dim, d_ocp_qp_sol& qp_sol, void * memory)
    {
        ::d_ocp_qp_sol_create(const_cast<d_ocp_qp_dim *>(&dim), &qp_sol, memory);
    }


    inline void ocp_qp_sol_get_all(d_ocp_qp_sol const& qp_sol, 
        double * const * u, double * const * x, double * const * ls, double * const * us,
        double * const * pi, double * const * lam_lb, double * const * lam_ub, double * const * lam_lg, double * const * lam_ug,
        double * const * lam_ls, double * const * lam_us)
    {
        ::d_ocp_qp_sol_get_all(const_cast<d_ocp_qp_sol *>(&qp_sol),
            const_cast<double **>(u), const_cast<double **>(x), const_cast<double **>(ls), const_cast<double **>(us),
            const_cast<double **>(pi), const_cast<double **>(lam_lb), const_cast<double **>(lam_ub), const_cast<double **>(lam_lg),
            const_cast<double **>(lam_ug), const_cast<double **>(lam_ls), const_cast<double **>(lam_us));
    }

    
    inline void ocp_qp_sol_get(char const * field, int stage, d_ocp_qp_sol const& qp_sol, double * vec)
    {
        ::d_ocp_qp_sol_get(const_cast<char *>(field), stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_get_u(int stage, d_ocp_qp_sol const& qp_sol, double *vec)
    {
        ::d_ocp_qp_sol_get_u(stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_get_x(int stage, d_ocp_qp_sol const& qp_sol, double *vec)
    {
        ::d_ocp_qp_sol_get_x(stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_get_sl(int stage, d_ocp_qp_sol const& qp_sol, double *vec)
    {
        ::d_ocp_qp_sol_get_sl(stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_get_su(int stage, d_ocp_qp_sol const& qp_sol, double *vec)
    {
        ::d_ocp_qp_sol_get_su(stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_get_pi(int stage, d_ocp_qp_sol const& qp_sol, double *vec)
    {
        ::d_ocp_qp_sol_get_pi(stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_get_lam_lb(int stage, d_ocp_qp_sol const& qp_sol, double *vec)
    {
        ::d_ocp_qp_sol_get_lam_lb(stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_get_lam_ub(int stage, d_ocp_qp_sol const& qp_sol, double *vec)
    {
        ::d_ocp_qp_sol_get_lam_ub(stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_get_lam_lg(int stage, d_ocp_qp_sol const& qp_sol, double *vec)
    {
        ::d_ocp_qp_sol_get_lam_lg(stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_get_lam_ug(int stage, d_ocp_qp_sol const& qp_sol, double *vec)
    {
        ::d_ocp_qp_sol_get_lam_ug(stage, const_cast<d_ocp_qp_sol *>(&qp_sol), vec);
    }

    
    inline void ocp_qp_sol_set(char const * field, int stage, double const * vec, d_ocp_qp_sol& qp_sol)
    {
        ::d_ocp_qp_sol_set(const_cast<char *>(field), stage, const_cast<double *>(vec), &qp_sol);
    }

    
    inline void ocp_qp_sol_set_u(int stage, double const * vec, d_ocp_qp_sol& qp_sol)
    {
        ::d_ocp_qp_sol_set_u(stage, const_cast<double *>(vec), &qp_sol);
    }

    
    inline void ocp_qp_sol_set_x(int stage, double const * vec, d_ocp_qp_sol& qp_sol)
    {
        ::d_ocp_qp_sol_set_x(stage, const_cast<double *>(vec), &qp_sol);
    }

    
    inline void ocp_qp_sol_set_sl(int stage, double const * vec, d_ocp_qp_sol& qp_sol)
    {
        ::d_ocp_qp_sol_set_sl(stage, const_cast<double *>(vec), &qp_sol);
    }

    
    inline void ocp_qp_sol_set_su(int stage, double const * vec, d_ocp_qp_sol& qp_sol)
    {
        ::d_ocp_qp_sol_set_su(stage, const_cast<double *>(vec), &qp_sol);
    }


    template <typename Real>
    struct OcpQpSol
    :   ocp_qp_sol<Real>
    {
        OcpQpSol(d_ocp_qp_dim const& dim)
        :   memory_ {new char[ocp_qp_sol_memsize(dim)]}
        {
            ocp_qp_sol_create(dim, *this, memory_.get());
        }


        void get_all(
            double * const * u, double * const * x, double * const * ls, double * const * us,
            double * const * pi, double * const * lam_lb, double * const * lam_ub, double * const * lam_lg, double * const * lam_ug,
            double * const * lam_ls, double * const * lam_us)
        {
            ocp_qp_sol_get_all(*this, 
                u, x, ls, us,
                pi, lam_lb, lam_ub, lam_lg, lam_ug,
                lam_ls, lam_us);
        }


        void get(char const * field, int stage, double * vec)
        {
            ocp_qp_sol_get(field, stage, *this, vec);
        }

        
        void get_u(int stage, double *vec) const
        {
            ocp_qp_sol_get_u(stage, *this, vec);
        }

        
        void get_x(int stage, double *vec) const
        {
            ocp_qp_sol_get_x(stage, *this, vec);
        }

        
        void get_sl(int stage, double *vec) const
        {
            ocp_qp_sol_get_sl(stage, *this, vec);
        }

        
        void get_su(int stage, double *vec) const
        {
            ocp_qp_sol_get_su(stage, *this, vec);
        }

        
        void get_pi(int stage, double *vec) const
        {
            ocp_qp_sol_get_pi(stage, *this, vec);
        }

        
        void get_lam_lb(int stage, double *vec) const
        {
            ocp_qp_sol_get_lam_lb(stage, *this, vec);
        }

        
        void get_lam_ub(int stage, double *vec) const
        {
            ocp_qp_sol_get_lam_ub(stage, *this, vec);
        }

        
        void get_lam_lg(int stage, double *vec) const
        {
            ocp_qp_sol_get_lam_lg(stage, *this, vec);
        }

        
        void get_lam_ug(int stage, double *vec) const
        {
            ocp_qp_sol_get_lam_ug(stage, *this, vec);
        }

        
        void set(char const * field, int stage, double const * vec)
        {
            ocp_qp_sol_set(const_cast<char *>(field), stage, const_cast<double *>(vec), *this);
        }

        
        void set_u(int stage, double const * vec)
        {
            ocp_qp_sol_set_u(stage, const_cast<double *>(vec), *this);
        }

        
        void set_x(int stage, double const * vec)
        {
            ocp_qp_sol_set_x(stage, const_cast<double *>(vec), *this);
        }

        
        void set_sl(int stage, double const * vec)
        {
            ocp_qp_sol_set_sl(stage, const_cast<double *>(vec), *this);
        }

        
        void set_su(int stage, double const * vec)
        {
            ocp_qp_sol_set_su(stage, const_cast<double *>(vec), *this);
        }

    
    private:
        std::unique_ptr<char []> memory_;
    };
}