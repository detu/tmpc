// **************************************************************************
//
// ocp_qp_dim API
//
// **************************************************************************

#pragma once

#include <hpipm_d_ocp_qp_dim.h>

#include <memory>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct OcpQpDimType
    {
        using Type = ::d_ocp_qp_dim;
    };


    template <typename Real>
    using ocp_qp_dim = typename OcpQpDimType<Real>::Type;
    
    
    template <typename Real>
    inline int ocp_qp_dim_memsize(int N);

    
    template <>
    inline int ocp_qp_dim_memsize<double>(int N)
    {
        return ::d_ocp_qp_dim_memsize(N);
    }


    inline void ocp_qp_dim_create(int N, d_ocp_qp_dim& qp_dim, void * memory)
    {
        ::d_ocp_qp_dim_create(N, &qp_dim, memory);
    }


    inline void ocp_qp_dim_set_all(
        int const * nx, int const * nu, int const * nbx, int const * nbu, int const * ng,
        int const * nsbx, int const * nsbu, int const * nsg, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_all(
            const_cast<int *>(nx), const_cast<int *>(nu), const_cast<int *>(nbx), const_cast<int *>(nbu), const_cast<int *>(ng),
            const_cast<int *>(nsbx), const_cast<int *>(nsbu), const_cast<int *>(nsg), &dim);
    }


    inline void ocp_qp_dim_set(char const * field, int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set(const_cast<char *>(field), stage, value, &dim);
    }


    inline void ocp_qp_dim_set_nx(int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_nx(stage, value, &dim);
    }


    
    inline void ocp_qp_dim_set_nu(int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_nu(stage, value, &dim);
    }
    
    
    inline void ocp_qp_dim_set_nbx(int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_nbx(stage, value, &dim);
    }
    
    
    inline void ocp_qp_dim_set_nbu(int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_nbu(stage, value, &dim);
    }
    
    
    inline void ocp_qp_dim_set_ng(int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_ng(stage, value, &dim);
    }
    
    
    inline void ocp_qp_dim_set_ns(int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_ns(stage, value, &dim);
    }
    
    
    inline void ocp_qp_dim_set_nsbx(int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_nsbx(stage, value, &dim);
    }
    
    inline void ocp_qp_dim_set_nsbu(int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_nsbu(stage, value, &dim);
    }
    
    
    inline void ocp_qp_dim_set_nsg(int stage, int value, d_ocp_qp_dim& dim)
    {
        ::d_ocp_qp_dim_set_nsg(stage, value, &dim);
    }
    
    
    inline int ocp_qp_dim_get(d_ocp_qp_dim const& dim, char const * field, int stage)
    {
        int value = -1;
        ::d_ocp_qp_dim_get(const_cast<d_ocp_qp_dim *>(&dim), const_cast<char *>(field), stage, &value);

        return value;
    }
    
    
    inline int ocp_qp_dim_get_N(d_ocp_qp_dim const& dim)
    {
        int value = -1;
        ::d_ocp_qp_dim_get_N(const_cast<d_ocp_qp_dim *>(&dim), &value);

        return value;
    }
    
    
    inline int ocp_qp_dim_get_nx(struct d_ocp_qp_dim *dim, int stage)
    {
        int value = -1;
        ::d_ocp_qp_dim_get_nx(const_cast<d_ocp_qp_dim *>(dim), stage, &value);

        return value;
    }
    
    
    inline int ocp_qp_dim_get_nu(struct d_ocp_qp_dim *dim, int stage)
    {
        int value = -1;
        ::d_ocp_qp_dim_get_nu(const_cast<d_ocp_qp_dim *>(dim), stage, &value);

        return value;
    }


    template <typename Real>
    struct OcpQpDim
    :	ocp_qp_dim<Real>
    {
        OcpQpDim(int N)
        :	memory_ {new char[ocp_qp_dim_memsize<Real>(N)]}
        {
            ocp_qp_dim_create(N, *this, memory_.get());
        }


        OcpQpDim(OcpQpDim&&) = default;


        void set_nx(int stage, int value)
        {
            ocp_qp_dim_set_nx(stage, value, *this);
        }


        void set_nu(int stage, int value)
        {
            ocp_qp_dim_set_nu(stage, value, *this);
        }


        /// @brief @brief Set number of NBX (state bound constraints)
        void set_nbx(int stage, int value)
        {
            ocp_qp_dim_set_nbx(stage, value, *this);
        }


        /// @brief Set number of NBU (input bound constraints)
        void set_nbu(int stage, int value)
        {
            ocp_qp_dim_set_nbu(stage, value, *this);
        }


        /// @brief Set number of NG (path constraints)
        void set_ng(int stage, int value)
        {
            ocp_qp_dim_set_ng(stage, value, *this);
        }


        void set_ns(int stage, int value)
        {
            ocp_qp_dim_set_ns(stage, value, *this);
        }


        /// @brief Set number of NSBX (soft state bound constraints)
        void set_nsbx(int stage, int value)
        {
            ocp_qp_dim_set_nsbx(stage, value, *this);
        }


        /// @brief Set number of NSBU (soft input bound constraints)
        void set_nsbu(int stage, int value)
        {
            ocp_qp_dim_set_nsbu(stage, value, *this);
        }


        /// @brief Set number of NSG (soft general bound constraints)
        void set_nsg(int stage, int value)
        {
            ocp_qp_dim_set_nsg(stage, value, *this);
        }


    private:
        std::unique_ptr<char []> memory_;
    };
}