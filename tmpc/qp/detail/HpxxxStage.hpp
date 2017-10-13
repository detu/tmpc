#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/ocp/OcpSolutionBase.hpp>
#include <tmpc/qp/OcpQpBase.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

#include "UnpaddedMatrix.hpp"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>


namespace tmpc :: detail
{
	///
	/// Common class for both HPMPC and HPIPM stage data.
	///
	template <typename Kernel_, StorageOrder SO>
	class HpxxxStage
	:	public OcpSolutionBase<HpxxxStage<Kernel_, SO>>
	,	public OcpQpBase<HpxxxStage<Kernel_, SO>>
	{
	public:
		using Kernel = Kernel_;
		using Real = typename Kernel::Real;

		HpxxxStage(OcpSize const& sz, size_t nx_next)
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
		,	A_ {nx_next, sz.nx()}
		,	B_ {nx_next, sz.nu()}
		,	b_(nx_next)
		,	C_ {sz.nc(), sz.nx()}
		,	D_ {sz.nc(), sz.nu()}
		,	lb_(sz.nu() + sz.nx())
		,	ub_(sz.nu() + sz.nx())
		,	lb_internal_(sz.nu() + sz.nx())
		,	ub_internal_(sz.nu() + sz.nx())
		,	lbd_(sz.nc())
		,	ubd_(sz.nc())
		,	idxs_(sz.ns())
		,	x_(sz.nx())
		,	u_(sz.nu())
		,	ls_(sz.ns())
		,	us_(sz.ns())
		,	pi_(nx_next)
		,	lam_(2 * sz.nc() + 2 * (sz.nx() + sz.nu()))
		,	lam_ls_(sz.ns())
		,	lam_us_(sz.ns())
		{
			// Initialize all numeric data to NaN so that if an uninitialized object
			// by mistake used in calculations is easier to detect.
			this->setNaN();

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

		HpxxxStage(HpxxxStage const&) = default;
		HpxxxStage(HpxxxStage &&) = default;

		auto const& A() const { return A_; }
		template <typename T> void A(T const& a) { noresize(A_) = a; }

		auto const& B() const { return B_; }
		template <typename T> void B(T const& b) { noresize(B_) = b; }

		auto const& b() const { return b_; }
		template <typename T> void b(T const& b) { noresize(b_) = b; }

		auto const& C() const { return C_; }
		template <typename T> void C(T const& c) { noresize(C_) = c; }

		auto const& D() const { return D_; }
		template <typename T> void D(T const& d) { noresize(D_) = d; }

		auto const& lbd() const {	return lbd_; }
		template <typename T> void lbd(T const& lbd) { noresize(lbd_) = lbd; }

		auto lbu() const { return subvector(lb_, 0, size_.nu());	}		
		
		// TODO: consider setting both upper and lower bounds at the same time.
		// Maybe create a Bounds class?
		template <typename T> void lbu(T const& lbu) 
		{ 
			subvector(lb_, 0, size_.nu()) = lbu;
		}

		auto lbx() const { return subvector(lb_, size_.nu(), size_.nx()); }
		template <typename T> void lbx(T const& lbx) 
		{ 
			subvector(lb_, size_.nu(), size_.nx()) = lbx; 
		}

		auto const& Q() const { return Q_; }
		template <typename T> void Q(T const& q) { noresize(Q_) = q; }

		auto const& R() const { return R_; }
		template <typename T> void R(T const& r) { noresize(R_) = r; }

		// HPMPC convention for S is [nu, nx], therefore the trans().
		auto S() const { return trans(S_); }
		void S(Real v) { S_ = v; }
		template <typename T> void S(T const& s) { noresize(S_) = trans(s); }

		auto const& q() const { return q_; }
		template <typename T> void q(T const& q) { noresize(q_) = q; }

		auto const& r() const { return r_; }
		template <typename T> void r(T const& r) { noresize(r_) = r; }

		// ----------------------------------------------
		// Soft constraints cost
		// ----------------------------------------------
		auto const& impl_Zl() const { return Zl_; }
		template <typename T> void impl_Zl(T const& val) { noresize(Zl_) = val; }

		auto const& impl_Zu() const { return Zu_; }
		template <typename T> void impl_Zu(T const& val) { noresize(Zu_) = val; }

		auto const& impl_zl() const { return zl_; }
		template <typename T> void impl_zl(T const& val) { noresize(zl_) = val; }

		auto const& impl_zu() const { return zu_; }
		template <typename T> void impl_zu(T const& val) { noresize(zu_) = val; }

		// ----------------------------------------------
		// Soft constraints index
		// ----------------------------------------------
		auto impl_idxs() const
		{
			return boost::make_iterator_range(idxs_);
		}


		template <typename IteratorRange>
		void impl_idxs(IteratorRange const& val)
		{
			if (val.size() != idxs_.size())
				throw std::invalid_argument("Soft constraints index size does not match");

			std::copy(val.begin(), val.end(), idxs_.begin());
		}

		// ----------------------------------------------
		auto const& ubd() const { return ubd_; }
		template <typename T> void ubd(T const& ubd) { noresize(ubd_) = ubd; }

		auto ubu() const { return subvector(ub_, 0, size_.nu()); }
		template <typename T> void ubu(T const& ubu) 
		{ 
			subvector(ub_, 0, size_.nu()) = ubu;
		}

		auto ubx() const { return subvector(ub_, size_.nu(), size_.nx()); }
		template <typename T> void ubx(T const& ubx) 
		{ 
			subvector(ub_, size_.nu(), size_.nx()) = ubx; 
		}

		auto const& x() const { return x_; }
		auto const& u() const { return u_;	}
		auto const& pi() const	{ return pi_; }
		
		auto lam_lbu() const 
		{ 
			return subvector(lam_, 2 * size_.nc(), size_.nu()); 
		}

		auto lam_ubu() const 
		{ 
			return subvector(lam_, 2 * size_.nc() + size_.nu(), size_.nu()); 
		}

		auto lam_lbx() const 
		{ 
			return subvector(lam_, 2 * size_.nc() + 2 * size_.nu(), size_.nx()); 
		}

		auto lam_ubx() const 
		{ 
			return subvector(lam_, 2 * size_.nc() + 2 * size_.nu() + size_.nx(), size_.nx()); 
		}

		auto lam_lbd() const 
		{ 
			return subvector(lam_, 0, size_.nc()); 
		}

		auto lam_ubd() const 
		{ 
			return subvector(lam_, size_.nc(), size_.nc()); 
		}

		OcpSize const& size() const { return size_; }

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
						throw std::invalid_argument("And invalid QP bound is found. For HPMPC, "
							"the bounds should be either both finite or [-inf, inf]");
				}
			}
		}

		// ******************************************************
		//                HPMPC raw data interface.
		//
		// The prefixes before _data() correspond to the names of
		// the argument to c_order_d_ip_ocp_hard_tv().
		// ******************************************************
		Real const * A_data () const { return A_.data(); }
		Real const * B_data () const { return B_.data(); }
		Real const * b_data () const { return b_.data(); }
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
		int nb() const { return idxb_.size(); }

		Real * x_data() { return x_.data(); }
		Real * u_data() { return u_.data(); }
		Real * ls_data() { return ls_.data(); }
		Real * us_data() { return us_.data(); }
		Real * pi_data() { return pi_.data(); }
		Real * lam_data() { return lam_.data(); }
		Real * lam_lb_data() { return lam_.data(); }
		Real * lam_ub_data() { return lam_lb_data() + size_.nx() + size_.nu(); }
		Real * lam_lg_data() { return lam_ub_data() + size_.nx() + size_.nu(); }
		Real * lam_ug_data() { return lam_lg_data() + size_.nc(); }
		Real * lam_ls_data() { return lam_ls_.data(); }
		Real * lam_us_data() { return lam_us_.data(); }

	private:
		OcpSize size_;

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

		// Inter-stage equalities x_{k+1} = A x_k + B u_k + c_k
		UnpaddedMatrix<Kernel, SO> A_;
		UnpaddedMatrix<Kernel, SO> B_;
		DynamicVector<Kernel> b_;

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

		// Lower and upper bound arrays for HPMPC,
		// containing finite values only.
		std::vector<Real> lb_internal_;
		std::vector<Real> ub_internal_;

		// Solution
		DynamicVector<Kernel> x_;
		DynamicVector<Kernel> u_;
		DynamicVector<Kernel> ls_;
		DynamicVector<Kernel> us_;
		DynamicVector<Kernel> pi_;
		DynamicVector<Kernel> lam_;
		DynamicVector<Kernel> lam_ls_;
		DynamicVector<Kernel> lam_us_;
	};
}
