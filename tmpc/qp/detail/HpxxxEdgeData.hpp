#pragma once

#include <tmpc/ocp/DynamicOcpSize.hpp>

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>

#include "UnpaddedMatrix.hpp"

#include <limits>


namespace tmpc :: detail
{
	///
	/// Common class for both HPMPC and HPIPM stage data.
	///
	template <typename Kernel_, StorageOrder SO>
	class HpxxxEdgeData
	{
	public:
		using Kernel = Kernel_;
		using Real = typename Kernel::Real;

		HpxxxEdgeData(DynamicOcpSize const& sz_u, DynamicOcpSize const& sz_v)
		:	A_ {sz_v.nx(), sz_u.nx()}
		,	B_ {sz_v.nx(), sz_u.nu()}
		,	b_(sz_v.nx())
		,	pi_(sz_v.nx())
		{
			// Initialize all numeric data to NaN so that if an uninitialized object
			// by mistake used in calculations is easier to detect.
			A_ = sNaN<Real>();
			B_ = sNaN<Real>();
			b_ = sNaN<Real>();
			pi_ = sNaN<Real>();
		}

		HpxxxEdgeData(HpxxxEdgeData const&) = default;
		HpxxxEdgeData(HpxxxEdgeData &&) = default;

		
		// ******************************************************
		// HPMPC/HPIPM raw data interface.
		// ******************************************************
		Real const * A_data () const { return A_.data(); }
		Real const * B_data () const { return B_.data(); }
		Real const * b_data () const { return b_.data(); }
		Real * pi_data() { return pi_.data(); }


	private:
		// Inter-stage equalities x_{k+1} = A x_k + B u_k + c_k
		UnpaddedMatrix<Kernel, SO> A_;
		UnpaddedMatrix<Kernel, SO> B_;
		DynamicVector<Kernel> b_;

		// Solution lagrange multipliers
		DynamicVector<Kernel> pi_;
	};
}
