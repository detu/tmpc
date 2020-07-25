#pragma once

#include "ErrorInfo.hpp"
#include "TreeOcpQpKkt.hpp"

#include <hpipm_d_tree_ocp_qp_ipm.h>
#include <hpipm_common.h>

#include <tmpc/qp/QpSolverException.hpp>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct TreeOcpQpIpmWorkspaceType;


    template <>
    struct TreeOcpQpIpmWorkspaceType<double>
    {
        using Type = ::d_tree_ocp_qp_ipm_workspace;
    };


    template <typename Real>
    using tree_ocp_qp_ipm_workspace = TreeOcpQpIpmWorkspaceType<Real>::Type;


    
    inline int memsize_tree_ocp_qp_ipm(d_tree_ocp_qp_dim const& dim, d_tree_ocp_qp_ipm_arg const& arg)
    {
        return ::d_memsize_tree_ocp_qp_ipm(const_cast<d_tree_ocp_qp_dim *>(&dim), const_cast<d_tree_ocp_qp_ipm_arg *>(&arg));
    }

    
    inline void create_tree_ocp_qp_ipm(d_tree_ocp_qp_dim const& dim, d_tree_ocp_qp_ipm_arg const& arg, d_tree_ocp_qp_ipm_workspace& ws, void * mem)
    {
        ::d_create_tree_ocp_qp_ipm(const_cast<d_tree_ocp_qp_dim *>(&dim), const_cast<d_tree_ocp_qp_ipm_arg *>(&arg), &ws, mem);
    }

    
    inline int get_tree_ocp_qp_ipm_iter(d_tree_ocp_qp_ipm_workspace const& ws)
    {
        return ::d_get_tree_ocp_qp_ipm_iter(const_cast<d_tree_ocp_qp_ipm_workspace *>(&ws));
    }

    
    inline double get_tree_ocp_qp_ipm_res_stat(d_tree_ocp_qp_ipm_workspace const& ws)
    {
        return ::d_get_tree_ocp_qp_ipm_res_stat(const_cast<d_tree_ocp_qp_ipm_workspace *>(&ws));
    }

    
    inline double get_tree_ocp_qp_ipm_res_eq(d_tree_ocp_qp_ipm_workspace const& ws)
    {
        return ::d_get_tree_ocp_qp_ipm_res_eq(const_cast<d_tree_ocp_qp_ipm_workspace *>(&ws));
    }

    
    inline double get_tree_ocp_qp_ipm_res_ineq(d_tree_ocp_qp_ipm_workspace const& ws)
    {
        return ::d_get_tree_ocp_qp_ipm_res_ineq(const_cast<d_tree_ocp_qp_ipm_workspace *>(&ws));
    }

    
    inline double get_tree_ocp_qp_ipm_res_comp(d_tree_ocp_qp_ipm_workspace const& ws)
    {
        return ::d_get_tree_ocp_qp_ipm_res_comp(const_cast<d_tree_ocp_qp_ipm_workspace *>(&ws));
    }

    
    inline double const * get_tree_ocp_qp_ipm_stat(d_tree_ocp_qp_ipm_workspace const& ws)
    {
        return ::d_get_tree_ocp_qp_ipm_stat(const_cast<d_tree_ocp_qp_ipm_workspace *>(&ws));
    }

    
    inline void solve_tree_ocp_qp_ipm(d_tree_ocp_qp const& qp,
        d_tree_ocp_qp_sol& qp_sol, d_tree_ocp_qp_ipm_arg const& arg, d_tree_ocp_qp_ipm_workspace& ws)
    {
        int const status = ::d_solve_tree_ocp_qp_ipm(
            const_cast<d_tree_ocp_qp *>(&qp), &qp_sol, const_cast<d_tree_ocp_qp_ipm_arg *>(&arg), &ws);
        
        if (status != SUCCESS)
            TMPC_THROW_EXCEPTION(QpSolverException {} << StatusErrorInfo {status} << boost::errinfo_api_function {"d_solve_tree_ocp_qp_ipm"});
    }


    template <typename Real>
    struct TreeOcpQpIpmWorkspace
    :   tree_ocp_qp_ipm_workspace<Real>
    {
        TreeOcpQpIpmWorkspace(d_tree_ocp_qp_dim const& dim, d_tree_ocp_qp_ipm_arg const& arg)
        :   memory_ {new char[memsize_tree_ocp_qp_ipm(dim, arg)]}
        {
            create_tree_ocp_qp_ipm(dim, arg, *this, memory_.get());
        }


        int get_iter() const
        {
            return get_tree_ocp_qp_ipm_iter(*this);
        }


        void solve(d_tree_ocp_qp const& qp, d_tree_ocp_qp_sol& qp_sol, d_tree_ocp_qp_ipm_arg const& arg)
        {
            solve_tree_ocp_qp_ipm(qp, qp_sol, arg, *this);
        }


        void fact_solve_kkt_unconstr(d_tree_ocp_qp const& qp, d_tree_ocp_qp_sol& qp_sol, d_tree_ocp_qp_ipm_arg const& arg)
        {
            fact_solve_kkt_unconstr_tree_ocp_qp(qp, qp_sol, arg, *this);
        }


    private:
        std::unique_ptr<char []> memory_;
    };
}