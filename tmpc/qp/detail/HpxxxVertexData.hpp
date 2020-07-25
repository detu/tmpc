#pragma once

#include <tmpc/ocp/DynamicOcpSize.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

#include "UnpaddedMatrix.hpp"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>


namespace tmpc :: detail
{
	///
	/// Common class for both HPMPC and HPIPM vertex data.
	///
	template <typename Kernel_, StorageOrder SO>
	class HpxxxVertexData
	{
	public:
		using Kernel = Kernel_;
		using Real = typename Kernel::Real;

		HpxxxVertexData(DynamicOcpSize const& sz)
		:	size_(sz)
		,	idxb_(sz.nu() + sz.nx())
		,	Q_ {sz.nx(), sz.nx()}
		,	R_ {sz.nu(), sz.nu()}
		,	S_ {sz.nu(), sz.nx()}	// <-- HPMPC convention for S is [nu, nx] (the corresponding cost term is u' * S_{hpmpc} * x)
		,	q_(sz.nx())
		,	r_(sz.nu())
		,	Zl_ {sz.ns(), sz.ns()}
		,	Zu_ {sz.ns(), sz.ns()}
		,	zl_(sz.ns())
		,	zu_(sz.ns())
		,	C_ {sz.nc(), sz.nx()}
		,	D_ {sz.nc(), sz.nu()}
		,	lb_(sz.nu() + sz.nx())
		,	ub_(sz.nu() + sz.nx())
		,	lb_internal_(sz.nu() + sz.nx())
		,	ub_internal_(sz.nu() + sz.nx())
		,	lbd_(sz.nc())
		,	ubd_(sz.nc())
		,	idxs_(sz.ns())
		,	lbls_(sz.ns())
		,	lbus_(sz.ns())
		,	x_(sz.nx())
		,	u_(sz.nu())
		,	ls_(sz.ns())
		,	us_(sz.ns())
		,	lam_(2 * sz.nc() + 2 * (sz.nx() + sz.nu()))
		,	lam_ls_(sz.ns())
		,	lam_us_(sz.ns())
		{
			// Initialize all numeric data to NaN so that if an uninitialized object
			// by mistake used in calculations is easier to detect.
			Q_ = sNaN<Real>();
			R_ = sNaN<Real>();
			S_ = sNaN<Real>();
			q_ = sNaN<Real>();
			r_ = sNaN<Real>();
			Zl_ = sNaN<Real>();
			Zu_ = sNaN<Real>();
			zl_ = sNaN<Real>();
			zu_ = sNaN<Real>();
			C_ = sNaN<Real>();
			D_ = sNaN<Real>();
			lb_ = sNaN<Real>();
			ub_ = sNaN<Real>();
			fill(lb_internal_.begin(), lb_internal_.end(), sNaN<Real>());
			fill(ub_internal_.begin(), ub_internal_.end(), sNaN<Real>());
			lbd_ = sNaN<Real>();
			ubd_ = sNaN<Real>();
			x_ = sNaN<Real>();
			u_ = sNaN<Real>();
			ls_ = sNaN<Real>();
			us_ = sNaN<Real>();
			lam_ = sNaN<Real>();
			lam_ls_ = sNaN<Real>();
			lam_us_ = sNaN<Real>();

			// Set lower bound for both slacks to 0.
			lbls_ = 0.;
			lbus_ = 0.;

			// idxb is initialized to its maximum size, s.t. nb == nx + nu.
			// This is necessary so that the solver workspace memory is calculated as its maximum when allocated.
			{
				int n = 0;
				std::generate(idxb_.begin(), idxb_.end(), [&n] { return n++; });
			}

			// Initialize idxs to 0,1,2, ... .
			{
				int n = 0;
				std::generate(idxs_.begin(), idxs_.end(), [&n] { return n++; });
			}
		}

		HpxxxVertexData(HpxxxVertexData const&) = default;
		HpxxxVertexData(HpxxxVertexData &&) = default;

		
		// Adjust hidxb so to account for infs in state and input bounds.
		void adjustBoundsIndex()
		{
			// this will not change the capacity and the data() pointers should stay the same.
			idxb_.clear();
			lb_internal_.clear();
			ub_internal_.clear();

			// Cycle through the bounds and check for infinities
			for (size_t i = 0; i < size_.nu() + size_.nx(); ++i)
			{
				if (std::isfinite(lb_[i]) && std::isfinite(ub_[i]))
				{
					// If both bounds are finite, add i to the bounds index,
					// and copy values to the lb_internal_ and ub_internal_.
					idxb_.push_back(i);
					lb_internal_.push_back(lb_[i]);
					ub_internal_.push_back(ub_[i]);
				}
				else 
				{
					// Otherwise, check that the values are [-inf, inf]
					if (!(lb_[i] == -inf<Real>() && ub_[i] == inf<Real>()))
						throw std::invalid_argument("An invalid QP bound is found. For HPMPC/HPIPM, "
							"the bounds should be either both finite or [-inf, inf]");
				}
			}
		}

		// ******************************************************
		// HPMPC/HPIPM raw data interface.
		// ******************************************************
		Real const * Q_data () const { return Q_.data(); }
		Real const * S_data () const { return S_.data(); }
		Real const * R_data () const { return R_.data(); }
		Real const * q_data () const { return q_.data(); }
		Real const * r_data () const { return r_.data(); }
		Real const * Zl_data() const { return Zl_.data(); }
		Real const * Zu_data() const { return Zu_.data(); }
		Real const * zl_data() const { return zl_.data(); }
		Real const * zu_data() const { return zu_.data(); }
		Real const * lb_data() const { return lb_internal_.data(); }
		Real const * ub_data() const { return ub_internal_.data(); }
		Real const * C_data () const { return C_.data(); }
		Real const * D_data () const { return D_.data(); }
		Real const * lg_data() const { return lbd_.data(); }
		Real const * ug_data() const { return ubd_.data(); }
		int const * hidxb_data() const { return idxb_.data(); }
		int const * idxs_data() const { return idxs_.data(); }
		Real const * lbls_data() const { return lbls_.data(); }
		Real const * lbus_data() const { return lbus_.data(); }


		// ******************************************************
		// Mutable pointers for accessing the data.
		// ******************************************************
		Real * Q_data() { return Q_.data(); }
		Real * S_data() { return S_.data(); }
		Real * R_data() { return R_.data(); }
		Real * q_data() { return q_.data(); }
		Real * r_data() { return r_.data(); }
		//Real * lx_data() { return ; }
		//Real * ux_data() { return ; }
		Real * Zl_data() { return Zl_.data(); }
		Real * Zu_data() { return Zu_.data(); }
		Real * zl_data() { return zl_.data(); }
		Real * zu_data() { return zu_.data(); }
		Real * C_data() { return C_.data(); }
		Real * D_data() { return D_.data(); }
		Real * lg_data() { return lg_.data(); }
		Real * ug_data() { return ug_.data(); }
		

		int nb() const { return idxb_.size(); }

		
		/// \brief Number of state bound constraints.
		int nbx() const
		{
			decltype(auto) lbx = this->lbx();
			decltype(auto) ubx = this->ubx();

			int count = 0;

			// Cycle through the bounds and check for infinities
			for (size_t i = 0; i < size_.nx(); ++i)
			{
				if (std::isfinite(lbx[i]) && std::isfinite(ubx[i]))
					++count;
			}

			return count;
		}


		/// \brief Number of input bound constraints.
		int nbu() const
		{
			decltype(auto) lbu = this->lbu();
			decltype(auto) ubu = this->ubu();

			int count = 0;

			// Cycle through the bounds and check for infinities
			for (size_t i = 0; i < size_.nu(); ++i)
			{
				if (std::isfinite(lbu[i]) && std::isfinite(ubu[i]))
					++count;
			}

			return count;
		}


		/// \brief Number of soft state bound constraints.
		int nsbx() const
		{
			int count = 0;

			// Cycle through the bounds and check for infinities
			for (size_t i = 0; i < size_.ns(); ++i)
			{
				if (size_.nu() <= idxs_[i] && idxs_[i] < size_.nu() + size_.nx())
					++count;
			}

			return count;
		}


		/// \brief Number of soft input bound constraints.
		int nsbu() const
		{
			int count = 0;

			// Cycle through the bounds and check for infinities
			for (size_t i = 0; i < size_.ns(); ++i)
			{
				if (idxs_[i] < size_.nu())
					++count;
			}

			return count;
		}


		/// \brief Number of soft general bound constraints.
		int nsg() const
		{
			int count = 0;

			// Cycle through the bounds and check for infinities
			for (size_t i = 0; i < size_.ns(); ++i)
			{
				if (idxs_[i] >= size_.nu() + size_.nx())
					++count;
			}

			return count;
		}


		Real * x_data() { return x_.data(); }
		Real * u_data() { return u_.data(); }
		Real * ls_data() { return ls_.data(); }
		Real * us_data() { return us_.data(); }
		Real * lam_data() { return lam_.data(); }
		Real * lam_lb_data() { return lam_.data(); }
		Real * lam_ub_data() { return lam_lb_data() + size_.nx() + size_.nu(); }
		Real * lam_lg_data() { return lam_ub_data() + size_.nx() + size_.nu(); }
		Real * lam_ug_data() { return lam_lg_data() + size_.nc(); }
		Real * lam_ls_data() { return lam_ls_.data(); }
		Real * lam_us_data() { return lam_us_.data(); }

	private:
		DynamicOcpSize size_;

		// Index of bound-constrained variables
		std::vector<int> idxb_;

		// Hessian = [R, S; S', Q]
		UnpaddedMatrix<Kernel, SO> Q_;
		UnpaddedMatrix<Kernel, SO> R_;
		UnpaddedMatrix<Kernel, SO> S_;	

		// Gradient = [r; q]
		DynamicVector<Kernel> r_;
		DynamicVector<Kernel> q_;

		// Soft constraints Hessian = [Zl, 0; 0, Zu]
		UnpaddedMatrix<Kernel, SO> Zl_;
		UnpaddedMatrix<Kernel, SO> Zu_;

		// Soft constraints gradient = [zl; zu]
		DynamicVector<Kernel> zl_;
		DynamicVector<Kernel> zu_;

		// Inequality constraints d_{min} <= C x_k + D u_k <= d_{max}
		UnpaddedMatrix<Kernel, SO> C_;
		UnpaddedMatrix<Kernel, SO> D_;
		DynamicVector<Kernel> lbd_;
		DynamicVector<Kernel> ubd_;

		// Bound constraints:
		// lb <= [u; x] <= ub
		DynamicVector<Kernel> lb_;
		DynamicVector<Kernel> ub_;

		// Soft constraints index
		std::vector<int> idxs_;

		// Lower bound of lower slack
		DynamicVector<Kernel> lbls_;

		// Lower bound of upper slack
		DynamicVector<Kernel> lbus_;

		// Lower and upper bound arrays for HPMPC/HPIPM,
		// containing finite values only.
		std::vector<Real> lb_internal_;
		std::vector<Real> ub_internal_;

		// Solution
		DynamicVector<Kernel> x_;
		DynamicVector<Kernel> u_;
		DynamicVector<Kernel> ls_;
		DynamicVector<Kernel> us_;
		DynamicVector<Kernel> lam_;
		DynamicVector<Kernel> lam_ls_;
		DynamicVector<Kernel> lam_us_;
	};
}
