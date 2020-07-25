#pragma once

#include <hpipm_d_tree_ocp_qp_kkt.h>


namespace tmpc :: hpipm
{
    inline void fact_solve_kkt_unconstr_tree_ocp_qp(d_tree_ocp_qp const& qp, 
        d_tree_ocp_qp_sol& qp_sol, d_tree_ocp_qp_ipm_arg const& arg, d_tree_ocp_qp_ipm_workspace& ws)
    {
        ::d_fact_solve_kkt_unconstr_tree_ocp_qp(const_cast<d_tree_ocp_qp *>(&qp), 
            &qp_sol, const_cast<d_tree_ocp_qp_ipm_arg *>(&arg), &ws);
    }
}