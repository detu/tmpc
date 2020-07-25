//  **************************************************************************
// 
//  ocp_qp API
// 
//  **************************************************************************
#pragma once

#include "OcpQpDim.hpp"

#include <hpipm_d_ocp_qp.h>

#include <memory>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct OcpQpType;


    template <>
    struct OcpQpType<double>
    {
        using Type = ::d_ocp_qp;
    };


    template <typename Real>
    using ocp_qp = typename OcpQpType<Real>::Type;

        
    inline int ocp_qp_memsize(d_ocp_qp_dim const& dim)
    {
        return ::d_ocp_qp_memsize(const_cast<d_ocp_qp_dim *>(&dim));
    }


    inline void ocp_qp_create(d_ocp_qp_dim const& dim, d_ocp_qp& qp, void * memory)
    {
        ::d_ocp_qp_create(const_cast<d_ocp_qp_dim *>(&dim), &qp, memory);
    }


    inline void ocp_qp_set_all(
        double const * const * A, double const * const * B, double const * const * b,
        double const * const * Q, double const * const * S, double const * const * R, double const * const * q, double const * const * r,
        int const * const * idxbx, double const * const * lbx, double const * const * ubx,
        int const * const * idxbu, double const * const * lbu, double const * const * ubu,
        double const * const * C, double const * const * D, double const * const * lg, double const * const * ug,
        double const * const * Zl, double const * const * Zu, double const * const * zl, double const * const * zu,
        int const * const * idxs, double const * const * ls, double const * const * us, d_ocp_qp& qp)
    {
        ::d_ocp_qp_set_all(
            const_cast<double **>(A), const_cast<double **>(B), const_cast<double **>(b), 
            const_cast<double **>(Q), const_cast<double **>(S), const_cast<double **>(R), const_cast<double **>(q), const_cast<double **>(r), 
            const_cast<int **>(idxbx), const_cast<double **>(lbx), const_cast<double **>(ubx),
            const_cast<int **>(idxbu), const_cast<double **>(lbu), const_cast<double **>(ubu),
            const_cast<double **>(C), const_cast<double **>(D), const_cast<double **>(lg), const_cast<double **>(ug),
            const_cast<double **>(Zl), const_cast<double **>(Zu), const_cast<double **>(zl), const_cast<double **>(zu),
            const_cast<int **>(idxs), const_cast<double **>(ls), const_cast<double **>(us), &qp);
    }


    inline void ocp_qp_set_A(int stage, double const * mat, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_A(stage, const_cast<double *>(mat), &qp);
    }
    
    
    inline void ocp_qp_set_B(int stage, double const * mat, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_B(stage, const_cast<double *>(mat), &qp);
    }
    
    
    inline void ocp_qp_set_b(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_b(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_Q(int stage, double const * mat, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_Q(stage, const_cast<double *>(mat), &qp);
    }
    
    
    inline void ocp_qp_set_S(int stage, double const * mat, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_S(stage, const_cast<double *>(mat), &qp);
    }
    
    
    inline void ocp_qp_set_R(int stage, double const * mat, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_R(stage, const_cast<double *>(mat), &qp);
    }
    
    
    inline void ocp_qp_set_q(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_q(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_r(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_r(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lb(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lb(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lb_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lb_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_ub(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_ub(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_ub_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_ub_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lbx(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lbx(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lbx_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lbx_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_el_lbx(int stage, int index, double elem, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_el_lbx(stage, index, &elem, &qp);
    }
    
    
    inline void ocp_qp_set_ubx(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_ubx(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_ubx_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_ubx_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_el_ubx(int stage, int index, double const * elem, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_el_ubx(stage, index, const_cast<double *>(elem), &qp);
    }
    
    
    inline void ocp_qp_set_lbu(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lbu(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lbu_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lbu_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_ubu(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_ubu(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_ubu_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_ubu_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_idxb(int stage, int const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_idxb(stage, const_cast<int *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_idxbx(int stage, int const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_idxbx(stage, const_cast<int *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_Jbx(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_Jbx(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_idxbu(int stage, int const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_idxbu(stage, const_cast<int *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_Jbu(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_Jbu(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_C(int stage, double const * mat, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_C(stage, const_cast<double *>(mat), &qp);
    }
    
    
    inline void ocp_qp_set_D(int stage, double const * mat, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_D(stage, const_cast<double *>(mat), &qp);
    }
    
    
    inline void ocp_qp_set_lg(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lg(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lg_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lg_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_ug(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_ug(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_ug_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_ug_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_Zl(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_Zl(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_Zu(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_Zu(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_zl(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_zl(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_zu(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_zu(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_idxs(int stage, int const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_idxs(stage, const_cast<int *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lls(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lls(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lls_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lls_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lus(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lus(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_lus_mask(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_lus_mask(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_Jsbu(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_Jsbu(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_Jsbx(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_Jsbx(stage, const_cast<double *>(vec), &qp);
    }
    
    
    inline void ocp_qp_set_Jsg(int stage, double const * vec, d_ocp_qp &qp)
    {
        ::d_ocp_qp_set_Jsg(stage, const_cast<double *>(vec), &qp);
    }
    

    template <typename Real>
    struct OcpQp
    :   ocp_qp<Real>
    {
        OcpQp(OcpQpDim<Real> const& dim)
        :   memory_ {new char[ocp_qp_memsize(dim)]}
        {
            ocp_qp_create(dim, *this, memory_.get());
        }

        OcpQp(OcpQp&&) = default;


        void set_A(int stage, double const * mat)
        {
            ocp_qp_set_A(stage, mat, *this);
        }
        
        
        void set_B(int stage, double const * mat)
        {
            ocp_qp_set_B(stage, mat, *this);
        }
        
        
        void set_b(int stage, double const * vec)
        {
            ocp_qp_set_b(stage, vec, *this);
        }
        
        
        void set_Q(int stage, double const * mat)
        {
            ocp_qp_set_Q(stage, mat, *this);
        }
        
        
        void set_S(int stage, double const * mat)
        {
            ocp_qp_set_S(stage, mat, *this);
        }
        
        
        void set_R(int stage, double const * mat)
        {
            ocp_qp_set_R(stage, mat, *this);
        }
        
        
        void set_q(int stage, double const * vec)
        {
            ocp_qp_set_q(stage, vec, *this);
        }
        
        
        void set_r(int stage, double const * vec)
        {
            ocp_qp_set_r(stage, vec, *this);
        }
        
        
        void set_lb(int stage, double const * vec)
        {
            ocp_qp_set_lb(stage, vec, *this);
        }
        
        
        void set_lb_mask(int stage, double const * vec)
        {
            ocp_qp_set_lb_mask(stage, vec, *this);
        }
        
        
        void set_ub(int stage, double const * vec)
        {
            ocp_qp_set_ub(stage, vec, *this);
        }
        
        
        void set_ub_mask(int stage, double const * vec)
        {
            ocp_qp_set_ub_mask(stage, vec, *this);
        }
        
        
        void set_lbx(int stage, double const * vec)
        {
            ocp_qp_set_lbx(stage, vec, *this);
        }
        
        
        void set_lbx_mask(int stage, double const * vec)
        {
            ocp_qp_set_lbx_mask(stage, vec, *this);
        }
        
        
        void set_el_lbx(int stage, int index, double const * elem)
        {
            ocp_qp_set_el_lbx(stage, elem *this);
        }
        
        
        void set_ubx(int stage, double const * vec)
        {
            ocp_qp_set_ubx(stage, vec, *this);
        }
        
        
        void set_ubx_mask(int stage, double const * vec)
        {
            ocp_qp_set_ubx_mask(stage, vec, *this);
        }
        
        
        void set_el_ubx(int stage, int index, double const * elem)
        {
            ocp_qp_set_el_ubx(stage, index, elem, *this);
        }
        
        
        void set_lbu(int stage, double const * vec)
        {
            ocp_qp_set_lbu(stage, vec, *this);
        }
        
        
        void set_lbu_mask(int stage, double const * vec)
        {
            ocp_qp_set_lbu_mask(stage, vec, *this);
        }
        
        
        void set_ubu(int stage, double const * vec)
        {
            ocp_qp_set_ubu(stage, vec, *this);
        }
        
        
        void set_ubu_mask(int stage, double const * vec)
        {
            ocp_qp_set_ubu_mask(stage, vec, *this);
        }
        
        
        void set_idxb(int stage, int const * vec)
        {
            ocp_qp_set_idxb(stage, vec, *this);
        }
        
        
        void set_idxbx(int stage, int const * vec)
        {
            ocp_qp_set_idxbx(stage, vec, *this);
        }
        
        
        void set_Jbx(int stage, double const * vec)
        {
            ocp_qp_set_Jbx(stage, vec, *this);
        }
        
        
        void set_idxbu(int stage, int const * vec)
        {
            ocp_qp_set_idxbu(stage, vec, *this);
        }
        
        
        void set_Jbu(int stage, double const * vec)
        {
            ocp_qp_set_Jbu(stage, vec, *this);
        }
        
        
        void set_C(int stage, double const * mat)
        {
            ocp_qp_set_C(stage, mat, *this);
        }
        
        
        void set_D(int stage, double const * mat)
        {
            ocp_qp_set_D(stage, mat, *this);
        }
        
        
        void set_lg(int stage, double const * vec)
        {
            ocp_qp_set_lg(stage, vec, *this);
        }
        
        
        void set_lg_mask(int stage, double const * vec)
        {
            ocp_qp_set_lg_mask(stage, vec, *this);
        }
        
        
        void set_ug(int stage, double const * vec)
        {
            ocp_qp_set_ug(stage, vec, *this);
        }
        
        
        void set_ug_mask(int stage, double const * vec)
        {
            ocp_qp_set_ug_mask(stage, vec, *this);
        }
        
        
        void set_Zl(int stage, double const * vec)
        {
            ocp_qp_set_Zl(stage, vec, *this);
        }
        
        
        void set_Zu(int stage, double const * vec)
        {
            ocp_qp_set_Zu(stage, vec, *this);
        }
        
        
        void set_zl(int stage, double const * vec)
        {
            ocp_qp_set_zl(stage, vec, *this);
        }
        
        
        void set_zu(int stage, double const * vec)
        {
            ocp_qp_set_zu(stage, vec, *this);
        }
        
        
        void set_idxs(int stage, int const * vec)
        {
            ocp_qp_set_idxs(stage, vec, *this);
        }
        
        
        void set_lls(int stage, double const * vec)
        {
            ocp_qp_set_lls(stage, vec, *this);
        }
        
        
        void set_lls_mask(int stage, double const * vec)
        {
            ocp_qp_set_lls_mask(stage, vec, *this);
        }
        
        
        void set_lus(int stage, double const * vec)
        {
            ocp_qp_set_lus(stage, vec, *this);
        }
        
        
        void set_lus_mask(int stage, double const * vec)
        {
            ocp_qp_set_lus_mask(stage, vec, *this);
        }
        
        
        void set_Jsbu(int stage, double const * vec)
        {
            ocp_qp_set_Jsbu(stage, vec, *this);
        }
        
        
        void set_Jsbx(int stage, double const * vec)
        {
            ocp_qp_set_Jsbx(stage, vec, *this);
        }
        
        
        void set_Jsg(int stage, double const * vec)
        {
            ocp_qp_set_Jsg(stage, vec, *this);
        }
    

    private:
        std::unique_ptr<char []> memory_;
    };
}