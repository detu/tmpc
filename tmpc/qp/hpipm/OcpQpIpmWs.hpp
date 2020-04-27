// **************************************************************************
//
// ocp_qp_ipm_ws API
//
// **************************************************************************
#pragma once

#include <tmpc/qp/QpSolverException.hpp>
#include <tmpc/qp/hpipm/ErrorInfo.hpp>

#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_common.h>

#include <memory>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct OcpQpIpmWsType;


    template <>
    struct OcpQpIpmWsType<double>
    {
        using Type = ::d_ocp_qp_ipm_ws;
    };


    template <typename Real>
    using ocp_qp_ipm_ws = typename OcpQpIpmWsType<Real>::Type;


    inline int ocp_qp_ipm_ws_memsize(d_ocp_qp_dim const& ocp_dim, d_ocp_qp_ipm_arg const& arg)
    {
        return ::d_ocp_qp_ipm_ws_memsize(const_cast<d_ocp_qp_dim *>(&ocp_dim), const_cast<d_ocp_qp_ipm_arg *>(&arg));
    }


    inline void ocp_qp_ipm_ws_create(d_ocp_qp_dim const& ocp_dim, d_ocp_qp_ipm_arg const& arg, d_ocp_qp_ipm_ws& ws, void * mem)
    {
        ::d_ocp_qp_ipm_ws_create(const_cast<d_ocp_qp_dim *>(&ocp_dim), const_cast<d_ocp_qp_ipm_arg *>(&arg), &ws, mem);
    }


    inline int ocp_qp_ipm_get_status(d_ocp_qp_ipm_ws const& ws)
    {
        int status = -1;
        ::d_ocp_qp_ipm_get_status(const_cast<d_ocp_qp_ipm_ws *>(&ws), &status);

        return status;
    }


    inline int ocp_qp_ipm_get_iter(d_ocp_qp_ipm_ws const& ws)
    {
        int iter = -1;
        ::d_ocp_qp_ipm_get_iter(const_cast<d_ocp_qp_ipm_ws *>(&ws), &iter);

        return iter;
    }


    inline void ocp_qp_ipm_solve(d_ocp_qp const& qp, d_ocp_qp_sol& qp_sol, d_ocp_qp_ipm_arg const& arg, d_ocp_qp_ipm_ws& ws)
    {
        ::d_ocp_qp_ipm_solve(const_cast<d_ocp_qp *>(&qp), &qp_sol, const_cast<d_ocp_qp_ipm_arg *>(&arg), &ws);
        
        int const status = ocp_qp_ipm_get_status(ws);
        if (status != SUCCESS)
            TMPC_THROW_EXCEPTION(QpSolverException {} << StatusErrorInfo {status} << boost::errinfo_api_function {"d_ocp_qp_ipm_solve"});
    }


    template <typename Real>
    struct OcpQpIpmWs
    :   ocp_qp_ipm_ws<Real>
    {
        OcpQpIpmWs(d_ocp_qp_dim const& ocp_dim, d_ocp_qp_ipm_arg const& arg)
        :   memory_ {new char[ocp_qp_ipm_ws_memsize(ocp_dim, arg)]}
        {
            ocp_qp_ipm_ws_create(ocp_dim, arg, *this, memory_.get());
        }


        int get_iter() const
        {
            return ocp_qp_ipm_get_iter(*this);
        }


    private:
        std::unique_ptr<char []> memory_;
    };
}