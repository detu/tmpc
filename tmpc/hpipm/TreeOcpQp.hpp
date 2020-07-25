#pragma once

#include <tmpc/hpipm/TreeOcpQpDim.hpp>

#include <hpipm_d_tree_ocp_qp.h>

#include <blasfeo/BlazeInterop.hpp>

#include <memory>
#include <cassert>


namespace tmpc :: hpipm
{
    template <typename Real>
    struct TreeOcpQpType;


    template <>
    struct TreeOcpQpType<double>
    {
        using Type = ::d_tree_ocp_qp;
    };


    template <typename Real>
    using tree_ocp_qp = TreeOcpQpType<Real>::Type;
    
    
    inline int memsize_tree_ocp_qp(d_tree_ocp_qp_dim const& dim)
    {
        return ::d_memsize_tree_ocp_qp(const_cast<d_tree_ocp_qp_dim *>(&dim));
    }
    
    
    inline void create_tree_ocp_qp(d_tree_ocp_qp_dim const& dim, d_tree_ocp_qp& qp, void * memory)
    {
        ::d_create_tree_ocp_qp(const_cast<d_tree_ocp_qp_dim *>(&dim), &qp, memory);
    }
    
    
    inline void cvt_colmaj_to_tree_ocp_qp(
        double const * const * A, double const * const * B, double const * const * b,
        double const * const * Q, double const * const * S, double const * const * R, double const * const * q, double const * const * r,
        int const * const * idxb, double const * const * d_lb, double const * const * d_ub,
        double const * const * C, double const * const * D, double const * const * d_lg, double const * const * d_ug,
        double const * const * Zl, double const * const * Zu, double const * const * zl, double const * const * zu,
        int const * const * idxs, double const * const * d_ls, double const * const * d_us, d_tree_ocp_qp& qp)
    {
        ::d_cvt_colmaj_to_tree_ocp_qp(
            const_cast<double **>(A), const_cast<double **>(B), const_cast<double **>(b),
            const_cast<double **>(Q), const_cast<double **>(S), const_cast<double **>(R), const_cast<double **>(q), const_cast<double **>(r),
            const_cast<int **>(idxb), const_cast<double **>(d_lb), const_cast<double **>(d_ub),
            const_cast<double **>(C), const_cast<double **>(D), const_cast<double **>(d_lg), const_cast<double **>(d_ug),
            const_cast<double **>(Zl), const_cast<double **>(Zu), const_cast<double **>(zl), const_cast<double **>(zu),
            const_cast<int **>(idxs), const_cast<double **>(d_ls), const_cast<double **>(d_us), &qp);
    }


    /// @brief Wrapper for tree_ocp_qp struct.
    ///
    /// The implementation of setter functions is based on cvt_colmaj_to_tree_ocp_qp() function,
    /// see https://github.com/giaf/hpipm/blob/master/tree_ocp_qp/x_tree_ocp_qp.c
    ///
    /// The layout of the bounds vector d is as follows:
    /// +---------+----+---------+----+
    /// |   lb    | lg |   ub    | ug |
    /// +---------+    +---------+    + 
    /// | lu | lx |    | uu | ux |    |
    /// +---------+----+---------+----+
    template <typename Real>
    struct TreeOcpQp
    :   tree_ocp_qp<Real>
    {
        /// @brief Construct tree_ocp_qp with given dimensions.
        ///
        /// @param dim Dimensions of the OCP QP. The lifetime of the dim object
        /// must exceed the lifetime of the TreeOcpQp being constructed,
        /// because tree_ocp_qp holds a pointer to the dim.
        ///
        explicit TreeOcpQp(TreeOcpQpDim<Real> const& dim)
        :   memory_ {new char[memsize_tree_ocp_qp(dim)]}
        {
            create_tree_ocp_qp(dim, *this, memory_.get());

            for (int ii = 0; ii < this->dim->Nn; ++ii)
            {
                // Default initialization of idxb
                for (int jj = 0; jj < this->dim->nb[ii]; ++jj)
                    this->idxb[ii][jj] = jj;

                // Default initialization of idxs
                for (int jj = 0; jj < this->dim->ns[ii]; ++jj)
                    this->idxs[ii][jj] = jj;

                // We are not ready to work with problems with not all components constrained.
                assert(this->dim->nb[ii] == this->dim->nbu[ii] + this->dim->nbx[ii]);
                assert(this->dim->nbu[ii] == this->dim->nu[ii]);
                assert(this->dim->nbx[ii] == this->dim->nx[ii]);
            }
        }


        template <typename MT, bool SO>
        void set_Q(int ii, blaze::Matrix<MT, SO> const& val)
        {
            assert(rows(val) == this->dim->nx[ii] && columns(val) == this->dim->nx[ii]);

            // CVT_MAT2STRMAT(nx[ii], nx[ii], Q[ii], nx[ii], qp->RSQrq+ii, nu[ii], nu[ii]);
		    blasfeo::pack(~val, this->RSQrq[ii], this->dim->nu[ii], this->dim->nu[ii]);		
        }


        template <typename MT, bool SO>
        void set_R(int ii, blaze::Matrix<MT, SO> const& val)
        {
            assert(rows(val) == this->dim->nu[ii] && columns(val) == this->dim->nu[ii]);

            // CVT_MAT2STRMAT(nu[ii], nu[ii], R[ii], nu[ii], qp->RSQrq+ii, 0, 0);
            blasfeo::pack(~val, this->RSQrq[ii], 0, 0);
        }


        template <typename MT, bool SO>
        void set_S(int ii, blaze::Matrix<MT, SO> const& val)
        {
            assert(rows(val) == this->dim->nu[ii] && columns(val) == this->dim->nx[ii]);

            // CVT_TRAN_MAT2STRMAT(nu[ii], nx[ii], S[ii], nu[ii], qp->RSQrq+ii, nu[ii], 0);
		    blasfeo::pack(trans(~val), this->RSQrq[ii], this->dim->nu[ii], 0);
        }

		
        template <typename VT>
        void set_r(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->nu[ii]);

            // CVT_TRAN_MAT2STRMAT(nu[ii], 1, r[ii], nu[ii], qp->RSQrq+ii, nu[ii]+nx[ii], 0);
            blasfeo::pack(trans(~val), this->RSQrq[ii], this->dim->nu[ii] + this->dim->nx[ii], 0);

            // CVT_VEC2STRVEC(nu[ii], r[ii], qp->rqz+ii, 0);
            blasfeo::pack(~val, this->rqz[ii], 0);
        }
		
		
        template <typename VT>
        void set_q(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->nx[ii]);

            // CVT_TRAN_MAT2STRMAT(nx[ii], 1, q[ii], nx[ii], qp->RSQrq+ii, nu[ii]+nx[ii], nu[ii]);
            blasfeo::pack(trans(~val), this->RSQrq[ii], this->dim->nu[ii] + this->dim->nx[ii], this->dim->nu[ii]);

            // CVT_VEC2STRVEC(nx[ii], q[ii], qp->rqz+ii, nu[ii]);
            blasfeo::pack(~val, this->rqz[ii], this->dim->nu[ii]);
        }


        template <typename MT, bool SO>
        void set_A(int ii, blaze::Matrix<MT, SO> const& val)
        {
            assert(rows(val) == this->dim->nx[ii] && columns(val) == this->dim->nx[idxdad(ii)]);

            // CVT_TRAN_MAT2STRMAT(nx[idx], nx[idxdad], A[ii], nx[idx], qp->BAbt+ii, nu[idxdad], 0);
            blasfeo::pack(trans(~val), this->BAbt[ii], this->dim->nu[idxdad(ii)], 0);
        }


        template <typename MT, bool SO>
        void set_B(int ii, blaze::Matrix<MT, SO> const& val)
        {
            assert(rows(val) == this->dim->nx[ii] && columns(val) == this->dim->nu[idxdad(ii)]);

            // CVT_TRAN_MAT2STRMAT(nx[idx], nu[idxdad], B[ii], nx[idx], qp->BAbt+ii, 0, 0);
            blasfeo::pack(trans(~val), this->BAbt[ii], 0, 0);
        }


        template <typename VT>
        void set_b(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->nx[ii]);

            // CVT_TRAN_MAT2STRMAT(nx[idx], 1, this->b[ii], nx[idx], qp->BAbt+ii, nu[idxdad]+nx[idxdad], 0);
            blasfeo::pack(trans(~val), this->BAbt[ii], this->dim->nu[idxdad(ii)] + this->dim->nx[idxdad(ii)], 0);

		    // CVT_VEC2STRVEC(nx[idx], this->b[ii], qp->b+ii, 0);
            blasfeo::pack(~val, this->b[ii], 0);
        }


        template <typename VT>
        void set_lb(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->nb[ii]);

            // CVT_VEC2STRVEC(nb[ii], d_lb[ii], qp->d+ii, 0);
            blasfeo::pack(~val, this->d[ii], 0);

            // VECSE_LIBSTR(nb[ii], 0.0, qp->m+ii, 0);
            blasfeo::vecse(this->dim->nb[ii], 0., this->m[ii], 0);
        }


        template <typename VT>
        void set_ub(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->nb[ii]);

            // CVT_VEC2STRVEC(nb[ii], d_ub[ii], qp->d+ii, nb[ii]+ng[ii]);
            blasfeo::pack(~val, this->d[ii], this->dim->nb[ii] + this->dim->ng[ii]);

            // VECSC_LIBSTR(nb[ii], -1.0, qp->d+ii, nb[ii]+ng[ii]);
            blasfeo::vecsc(this->dim->nb[ii], -1., this->d[ii], this->dim->nb[ii] + this->dim->ng[ii]);

            // VECSE_LIBSTR(nb[ii], 0.0, qp->m+ii, nb[ii]+ng[ii]);
            blasfeo::vecse(this->dim->nb[ii], 0., this->m[ii], this->dim->nb[ii] + this->dim->ng[ii]);
        }


        template <typename VT>
        void set_lbu(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->nbu[ii]);
            setLowerBound(ii, val, 0);
        }


        template <typename VT>
        void set_ubu(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->nbu[ii]);
            setUpperBound(ii, val, 0);
        }


        template <typename VT>
        void set_lbx(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->nbx[ii]);
            setLowerBound(ii, val, this->dim->nbu[ii]);
        }


        template <typename VT>
        void set_ubx(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->nbx[ii]);
            setUpperBound(ii, val, this->dim->nbu[ii]);
        }


        template <typename MT, bool SO>
        void set_C(int ii, blaze::Matrix<MT, SO> const& val)
        {
            assert(rows(val) == this->dim->ng[ii] && columns(val) == this->dim->nx[ii]);

            // CVT_TRAN_MAT2STRMAT(ng[ii], nx[ii], C[ii], ng[ii], qp->DCt+ii, nu[ii], 0);
            blasfeo::pack(trans(~val), this->DCt[ii], this->dim->nu[ii], 0);
        }


        template <typename MT, bool SO>
        void set_D(int ii, blaze::Matrix<MT, SO> const& val)
        {
            assert(rows(val) == this->dim->ng[ii] && columns(val) == this->dim->nu[ii]);

            // CVT_TRAN_MAT2STRMAT(ng[ii], nu[ii], D[ii], ng[ii], qp->DCt+ii, 0, 0);
            blasfeo::pack(trans(~val), this->DCt[ii], 0, 0);
        }


        template <typename VT>
        void set_lg(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->ng[ii]);
            
            // CVT_VEC2STRVEC(ng[ii], d_lg[ii], qp->d+ii, nb[ii]);
            blasfeo::pack(~val, this->d[ii], this->dim->nb[ii]);

            // VECSE_LIBSTR(ng[ii], 0.0, qp->m+ii, nb[ii]);
            blasfeo::vecse(this->dim->ng[ii], 0., this->m[ii], this->dim->nb[ii]);
        }


        template <typename VT>
        void set_ug(int ii, blaze::Vector<VT, blaze::columnVector> const& val)
        {
            assert(size(val) == this->dim->ng[ii]);
            
            // CVT_VEC2STRVEC(ng[ii], d_ug[ii], qp->d+ii, 2*nb[ii]+ng[ii]);
            blasfeo::pack(~val, this->d[ii], 2 * this->dim->nb[ii] + this->dim->ng[ii]);

			// VECSC_LIBSTR(ng[ii], -1.0, qp->d+ii, 2*nb[ii]+ng[ii]);
            blasfeo::vecsc(this->dim->ng[ii], -1., this->d[ii], 2 * this->dim->nb[ii] + this->dim->ng[ii]);

            // VECSE_LIBSTR(ng[ii], 0.0, qp->m+ii, 2*nb[ii]+ng[ii]);
            blasfeo::vecse(this->dim->ng[ii], 0., this->m[ii], 2 * this->dim->nb[ii] + this->dim->ng[ii]);
        }


        void set_idxb(int ii, int const * idxb)
        {
            for (int jj = 0; jj < this->dim->nb[ii]; ++jj)
				this->idxb[ii][jj] = idxb[jj];
        }


        void set_idxbx(int stage, int const * idxbx)
        {
            // Implementation based on OCP_QP_SET_IDXBX(),
            // see https://github.com/giaf/hpipm/blob/master/ocp_qp/x_ocp_qp.c
            for (int ii = 0; ii < this->dim->nbx[stage]; ++ii)
                this->idxb[stage][this->dim->nbu[stage] + ii] = this->dim->nu[stage] + idxbx[ii];
        }


        void set_idxbu(int stage, int const * idxbu)
        {
            // Implementation based on OCP_QP_SET_IDXBU(),
            // see https://github.com/giaf/hpipm/blob/master/ocp_qp/x_ocp_qp.c
            for (int ii = 0; ii < this->dim->nbu[stage]; ++ii)
                this->idxb[stage][ii] = idxbu[ii];
        }


    private:
        std::unique_ptr<char []> memory_;


        auto idxdad(int ii) const noexcept
        {
            // idx = ii+1;
		    // idxdad = (ttree->root+idx)->dad;
            return this->dim->ttree->root[ii + 1].dad;
        }


        template <typename VT>
        void setLowerBound(int ii, blaze::Vector<VT, blaze::columnVector> const& val, int idx)
        {
            blasfeo::pack(~val, this->d[ii], idx);
            blasfeo::vecse(size(val), 0., this->m[ii], idx);
        }


        template <typename VT>
        void setUpperBound(int ii, blaze::Vector<VT, blaze::columnVector> const& val, int idx)
        {
            blasfeo::pack(~val, this->d[ii], this->dim->nb[ii] + this->dim->ng[ii] + idx);
            blasfeo::vecsc(size(val), -1., this->d[ii], this->dim->nb[ii] + this->dim->ng[ii] + idx);
            blasfeo::vecse(size(val), 0., this->m[ii], this->dim->nb[ii] + this->dim->ng[ii] + idx);
        }
    };
}