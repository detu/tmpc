// **************************************************************************
//
// ocp_qp_ipm_arg API
//
// **************************************************************************
#pragma once

#include <tmpc/qp/hpipm/OcpQpDim.hpp>

#include <hpipm_d_ocp_qp_ipm.h>

#include <memory>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct OcpQpIpmArgType;


    template <>
    struct OcpQpIpmArgType<double>
    {
        using Type = ::d_ocp_qp_ipm_arg;
    };


    template <typename Real>
    using ocp_qp_ipm_arg = typename OcpQpIpmArgType<Real>::Type;


    inline int ocp_qp_ipm_arg_memsize(d_ocp_qp_dim const& ocp_dim)
    {
        return ::d_ocp_qp_ipm_arg_memsize(const_cast<d_ocp_qp_dim *>(&ocp_dim));
    }


    inline void ocp_qp_ipm_arg_create(d_ocp_qp_dim const& ocp_dim, d_ocp_qp_ipm_arg& arg, void * mem)
    {
        ::d_ocp_qp_ipm_arg_create(const_cast<d_ocp_qp_dim *>(&ocp_dim), &arg, mem);
    }


    inline void ocp_qp_ipm_arg_set_default(hpipm_mode mode, d_ocp_qp_ipm_arg& arg)
    {
        ::d_ocp_qp_ipm_arg_set_default(mode, &arg);
    }


    /// @brief Set maximum number of iterations
    inline void ocp_qp_ipm_arg_set_iter_max(int iter_max, d_ocp_qp_ipm_arg& arg)
    {
        ::d_ocp_qp_ipm_arg_set_iter_max(&iter_max, &arg);
    }


    // @brief Set minimum step lenght
    inline void ocp_qp_ipm_arg_set_alpha_min(double alpha_min, d_ocp_qp_ipm_arg& arg)
    {
        ::d_ocp_qp_ipm_arg_set_alpha_min(&alpha_min, &arg);
    }
    
    
    template <typename Real>
    struct OcpQpIpmArg
    :	ocp_qp_ipm_arg<Real>
    {
        OcpQpIpmArg(ocp_qp_dim<Real> const& dim, hpipm_mode mode)
        :	memory_ {new char[ocp_qp_ipm_arg_memsize(dim)]}
        {
            ocp_qp_ipm_arg_create(dim, *this, memory_.get());
            ocp_qp_ipm_arg_set_default(mode, *this);
        }


        /// @brief Set maximum number of iterations
        void set_iter_max(int iter_max)
        {
            ocp_qp_ipm_arg_set_iter_max(iter_max, *this);
        }


        // @brief Set minimum step lenght
        void set_alpha_min(double alpha_min, d_ocp_qp_ipm_arg& arg)
        {
            ocp_qp_ipm_arg_set_alpha_min(alpha_min, *this);
        }


    private:
        std::unique_ptr<char []> memory_;
    };
}