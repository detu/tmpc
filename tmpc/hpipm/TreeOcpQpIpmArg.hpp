#pragma once

#include <hpipm_d_tree_ocp_qp_ipm.h>

#include <memory>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct TreeOcpQpIpmArgType;


    template <>
    struct TreeOcpQpIpmArgType<double>
    {
        using Type = ::d_tree_ocp_qp_ipm_arg;
    };


    template <typename Real>
    using tree_ocp_qp_ipm_arg = typename TreeOcpQpIpmArgType<Real>::Type;


    
    inline int memsize_tree_ocp_qp_ipm_arg(d_tree_ocp_qp_dim const& dim)
    {
        return ::d_memsize_tree_ocp_qp_ipm_arg(const_cast<d_tree_ocp_qp_dim *>(&dim));
    }

    
    inline void create_tree_ocp_qp_ipm_arg(d_tree_ocp_qp_dim const& dim, d_tree_ocp_qp_ipm_arg& arg, void * mem)
    {
        ::d_create_tree_ocp_qp_ipm_arg(const_cast<d_tree_ocp_qp_dim *>(&dim), &arg, mem);
    }


    inline void set_default_tree_ocp_qp_ipm_arg(hpipm_mode mode, d_tree_ocp_qp_ipm_arg& arg)
    {
        ::d_set_default_tree_ocp_qp_ipm_arg(mode, &arg);
    }

    
    inline void set_tree_ocp_qp_ipm_arg_iter_max(int iter_max, d_tree_ocp_qp_ipm_arg& arg)
    {
        ::d_set_tree_ocp_qp_ipm_arg_iter_max(iter_max, &arg);
    }

    
    inline void set_tree_ocp_qp_ipm_arg_mu0(double mu0, d_tree_ocp_qp_ipm_arg& arg)
    {
        ::d_set_tree_ocp_qp_ipm_arg_mu0(mu0, &arg);
    }

    
    inline void set_tree_ocp_qp_ipm_arg_tol_stat(double tol_stat, d_tree_ocp_qp_ipm_arg& arg)
    {
        ::d_set_tree_ocp_qp_ipm_arg_tol_stat(tol_stat, &arg);
    }

    
    inline void set_tree_ocp_qp_ipm_arg_tol_eq(double tol_eq, d_tree_ocp_qp_ipm_arg& arg)
    {
        ::d_set_tree_ocp_qp_ipm_arg_tol_eq(tol_eq, &arg);
    }
    

    inline void set_tree_ocp_qp_ipm_arg_tol_ineq(double tol_ineq, d_tree_ocp_qp_ipm_arg& arg)
    {
        ::d_set_tree_ocp_qp_ipm_arg_tol_ineq(tol_ineq, &arg);
    }

    
    inline void set_tree_ocp_qp_ipm_arg_tol_comp(double tol_comp, d_tree_ocp_qp_ipm_arg& arg)
    {
        ::d_set_tree_ocp_qp_ipm_arg_tol_comp(tol_comp, &arg);
    }

    
    inline void set_tree_ocp_qp_ipm_arg_reg_prim(double reg, d_tree_ocp_qp_ipm_arg& arg)
    {
        ::d_set_tree_ocp_qp_ipm_arg_reg_prim(reg, &arg);
    }


    template <typename Real>
    struct TreeOcpQpIpmArg
    :   tree_ocp_qp_ipm_arg<Real>
    {
        TreeOcpQpIpmArg(d_tree_ocp_qp_dim const& dim, hpipm_mode mode)
        :	memory_ {new char[memsize_tree_ocp_qp_ipm_arg(dim)]}
        {
            create_tree_ocp_qp_ipm_arg(dim, *this, memory_.get());
            set_default_tree_ocp_qp_ipm_arg(mode, *this);            
        }


        void set_iter_max(int iter_max)
        {
            set_tree_ocp_qp_ipm_arg_iter_max(iter_max, *this);
        }

        
        void set_mu0(double mu0)
        {
            set_tree_ocp_qp_ipm_arg_mu0(mu0, *this);
        }

        
        void set_tol_stat(double tol_stat)
        {
            set_tree_ocp_qp_ipm_arg_tol_stat(tol_stat, *this);
        }

        
        void set_tol_eq(double tol_eq)
        {
            set_tree_ocp_qp_ipm_arg_tol_eq(tol_eq, *this);
        }

        
        void set_tol_ineq(double tol_ineq)
        {
            set_tree_ocp_qp_ipm_arg_tol_ineq(tol_ineq, *this);
        }

        
        void set_tol_comp(double tol_comp)
        {
            set_tree_ocp_qp_ipm_arg_tol_comp(tol_comp, *this);
        }

        
        void set_reg_prim(double reg)
        {
            set_tree_ocp_qp_ipm_arg_reg_prim(reg, *this);
        }


    private:
        std::unique_ptr<char []> memory_;
    };
}