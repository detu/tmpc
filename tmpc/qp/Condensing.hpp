#pragma once

#include <tmpc/ocp/OcpSize.hpp>
#include <tmpc/qp/OcpQpBase.hpp>
#include <tmpc/Matrix.hpp>

#include <vector>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <array>

namespace tmpc
{
	/**
	 * Resulting QP size after full condensing.
	 */
	OcpSize condensedQpSize(std::vector<OcpSize> const& sz);

	template <typename InIter>
	OcpSize condensedQpSize(InIter sz_begin, InIter sz_end)
	{
		if (sz_begin == sz_end)
			throw std::invalid_argument("condensedQpSize(): OcpSize range must be not empty");

		return OcpSize(
			sz_begin->nx(),
			std::accumulate(sz_begin, sz_end, std::size_t{0},
				[] (std::size_t n, OcpSize const& s) { return n + s.nu(); }),
			std::accumulate(sz_begin, sz_end, std::size_t{0},
				[] (std::size_t n, OcpSize const& s) { return n + s.nc(); })
			+ std::accumulate(sz_begin + 1, sz_end, std::size_t{0},
					[] (std::size_t n, OcpSize const& s) { return n + s.nx(); })
			);
	}

	/**
	 * \brief Condensing algorithm.
	 *
	 * Manages resources needed for the condensing algorithm of a given size.
	 *
	 * \tparam <Kernel> Kernel for matrix arithmetic.
	 * 
	 * TODO: account for soft constraints.
	 */
	template <typename Kernel_>
	class Condensing
	{
	public:
		using Kernel = Kernel_;

		class CondensedStage
		:	public OcpQpBase<CondensedStage>
		{
			using size_type = typename Kernel::size_t;
			using Matrix = DynamicMatrix<Kernel>;
			using Vector = DynamicVector<Kernel>;

		public:
			CondensedStage(CondensedStage const&) = delete;
			CondensedStage(CondensedStage &&) = default;

			CondensedStage(Condensing const& c)
			:	c_(c)
			{
			}
	
			auto const& A() const {
				return c_.A_;
			}
	
			auto const& B() const {
				return c_.B_;
			}
	
			auto const& b() const {
				return c_.b_;
			}
	
			auto const& C() const {
				return c_.Cc_;
			}
	
			auto const& D() const {
				return c_.Dc_;
			}
		
			auto const& lbd() const {
				return c_.lbd_;
			}
	
			auto const& lbu() const {
				return c_.lbu_;
			}
	
			auto const& lbx() const {
				return c_.lbx_;
			}
	
			auto const& Q() const {
				return c_.Qc_;
			}
	
			auto const& R() const {
				return c_.Rc_;
			}
	
			auto const& S() const {
				return c_.Sc_;
			}
		
			auto const& q() const {
				return c_.qc_;
			}
	
			auto const& r() const {
				return c_.rc_;
			}
	
			auto const& ubd() const {
				return c_.ubd_;
			}
		
			auto const& ubu() const {
				return c_.ubu_;
			}
		
			auto const& ubx() const {
				return c_.ubx_;
			}


			// -----------------------------------------------------------
			// Soft constraints cost
			// -----------------------------------------------------------
			decltype(auto) impl_Zl() const 
			{ 
				return StaticMatrix<Kernel, 0, 0> {};
			}

			
			decltype(auto) impl_Zu() const 
			{	
				return StaticMatrix<Kernel, 0, 0> {};
			}

			
			decltype(auto) impl_zl() const 
			{ 
				return StaticVector<Kernel, 0> {}; 
			}

			
			decltype(auto) impl_zu() const 
			{	
				return StaticVector<Kernel, 0> {};
			}


			// -----------------------------------------------------------
			// Soft constraints index
			// -----------------------------------------------------------
			decltype(auto) impl_idxs() const
			{
				return std::array<int, 0> {};
			}


			// -----------------------------------------------------------	
			OcpSize const& size() const
			{
				return c_.cs_;
			}

			size_t nxNext() const
			{
				return rows(c_.A_);
			}
	
		private:
			Condensing const& c_;
		};

		template <typename InIter>
		Condensing(InIter const& sz_first, InIter const& sz_last)
		:	cs_(condensedQpSize(sz_first, sz_last))
		,	Qc_(cs_.nx(), cs_.nx())
		,   Rc_(cs_.nu(), cs_.nu())
		,   Sc_(cs_.nx(), cs_.nu())
		,   qc_(cs_.nx())
		,   rc_(cs_.nu())
		,   Cc_(cs_.nc(), cs_.nx())
		,   Dc_(cs_.nc(), cs_.nu(), 0)
		,   lbd_(cs_.nc())
		,   ubd_(cs_.nc())
		,   lbu_(cs_.nu())
		,   ubu_(cs_.nu())
		{
		}

		/**
		* \brief Perform the condensing.
		*/
		template <typename InIter>
		CondensedStage operator()(InIter first, InIter last)
		{
			if (condensedQpSize(qpSizeIterator(first), qpSizeIterator(last)) != cs_)
				throw std::invalid_argument("QP size does not match the condensed size");

			// Initialization of QP variables
			auto nu = first->size().nu();
			auto nc = first->size().nc();

			Qc_ = first->Q();
			submatrix(Rc_, 0, 0, nu, nu) = first->R();
			//K::left_cols(Sc_, nu) = first_->S();
			submatrix(Sc_, 0, 0, Sc_.rows(), nu) = first->S();
			qc_ = first->q();
			//K::head(rc_, nu) = first_->r();
			subvector(rc_, 0, nu) = first->r();

			//K::top_rows(Cc_, nc) = first_->C();
			submatrix(Cc_, 0, 0, nc, columns(Cc_)) = first->C();
			//K::top_left_corner(Dc_, nc, nu) = first_->D();
			submatrix(Dc_, 0, 0, nc, nu) = first->D();
			//K::head(lbd, nc) = first_->lbd();
			subvector(lbd_, 0, nc) = first->lbd();
			//K::head(ubd, nc) = first_->ubd();
			subvector(ubd_, 0, nc) = first->ubd();

			//K::head(lbu, nu) = first_->lbu();
			subvector(lbu_, 0, nu) = first->lbu();
			//K::head(ubu, nu) = first_->ubu();
			subvector(ubu_, 0, nu) = first->ubu();

			A_ = first->A();
			B_ = first->B();
			b_ = first->b();

			for (auto stage = first + 1; stage != last; ++stage)
			{
				auto const& sz = stage->size();
				// TODO: add nx1() to OcpSize
				auto const nx_next = stage->B().rows();

				// Update Q
				Qc_ += trans(A_) * stage->Q() * A_;

				// Update S
				//K::left_cols(Sc_, nu) += K::trans(A) * stage->Q() * B;
				submatrix(Sc_, 0, 0, Sc_.rows(), nu) += trans(A_) * stage->Q() * B_;
				//K::middle_cols(Sc_, nu, sz.nu()) = K::trans(A) * stage->S();
				submatrix(Sc_, 0, nu, Sc_.rows(), sz.nu()) = trans(A_) * stage->S();

				// Update R
				//K::top_left_corner(Rc_, nu, nu) += K::trans(B) * stage->Q() * B;
				submatrix(Rc_, 0, 0, nu, nu) += trans(B_) * stage->Q() * B_;
				//K::block(Rc_, 0, nu, nu, sz.nu()) = K::trans(B) * stage->S();
				submatrix(Rc_, 0, nu, nu, sz.nu()) = trans(B_) * stage->S();
				//K::block(Rc_, nu, 0, sz.nu(), nu) = K::trans(stage->S()) * B;
				submatrix(Rc_, nu, 0, sz.nu(), nu) = trans(stage->S()) * B_;
				//K::block(Rc_, nu, nu, sz.nu(), sz.nu()) = stage->R();
				submatrix(Rc_, nu, nu, sz.nu(), sz.nu()) = stage->R();

				// Update q
				qc_ += trans(A_) * (stage->Q() * b_ + stage->q());

				// Update r
				//K::head(rc_, nu) += K::trans(B) * (stage->Q() * b + stage->q());
				subvector(rc_, 0, nu) += trans(B_) * (stage->Q() * b_ + stage->q());
				//K::segment(rc_, nu, sz.nu()) = K::trans(stage->S()) * b + stage->r();
				subvector(rc_, nu, sz.nu()) = trans(stage->S()) * b_ + stage->r();

				// Update C
				//K::middle_rows(Cc_, nc, sz.nx() + sz.nc()) << A, stage->C() * A;
				submatrix(Cc_, nc          , 0, sz.nx(), columns(Cc_)) = A_;
				submatrix(Cc_, nc + sz.nx(), 0, sz.nc(), columns(Cc_)) = stage->C() * A_;

				// Update D (which is initialized to 0)
				//K::middle_rows(Dc_, nc, sz.nx() + sz.nc()) << B, K::zero(sz.nx(), cs_.nu() - nu),
				//									  stage->C() * B, stage->D(), K::zero(sz.nc(), cs_.nu() - nu - sz.nu());
				submatrix(Dc_, nc          ,  0, sz.nx(),      nu) = B_;
				submatrix(Dc_, nc + sz.nx(),  0, sz.nc(),      nu) = stage->C() * B_;
				submatrix(Dc_, nc + sz.nx(), nu, sz.nc(), sz.nu()) = stage->D();

				// Update lbd
				//K::segment(lbd, nc, sz.nx() + sz.nc()) << stage->lbx() - b, stage->lbd() - stage->C() * b;
				subvector(lbd_, nc          , sz.nx()) = stage->lbx() - b_;
				subvector(lbd_, nc + sz.nx(), sz.nc()) = stage->lbd() - stage->C() * b_;

				// Update ubd
				//K::segment(ubd, nc, sz.nx() + sz.nc()) << stage->ubx() - b, stage->ubd() - stage->C() * b;
				subvector(ubd_, nc          , sz.nx()) = stage->ubx() - b_;
				subvector(ubd_, nc + sz.nx(), sz.nc()) = stage->ubd() - stage->C() * b_;

				// Uplate lbu, ubu
				subvector(lbu_, nu, sz.nu()) = stage->lbu();
				subvector(ubu_, nu, sz.nu()) = stage->ubu();

				// Update A
				A_ = stage->A() * A_;

				// Update B
				{
					DynamicMatrix<Kernel> B_next(nx_next, nu + sz.nu());
					submatrix(B_next, 0,  0, nx_next,      nu) = stage->A() * B_;
					submatrix(B_next, 0, nu, nx_next, sz.nu()) = stage->B();
					B_ = std::move(B_next);
				}

				// Update b
				b_ = stage->A() * b_ + stage->b();

				// Update indices
				nu += sz.nu();
				nc += sz.nx() + sz.nc();
			}

			// Set the state bounds of the condensed stage
			lbx_ = first->lbx();
			ubx_ = first->ubx();

			CondensedStage result(*this);
			return std::move(result);
		}


		auto const& condensedSize() const
		{
			return cs_;
		}
		

	private:
		OcpSize const cs_;

		DynamicMatrix<Kernel> Qc_;
		DynamicMatrix<Kernel> Rc_;
		DynamicMatrix<Kernel> Sc_;
		DynamicVector<Kernel> qc_;
		DynamicVector<Kernel> rc_;

		DynamicMatrix<Kernel> Cc_;
		DynamicMatrix<Kernel> Dc_;
		DynamicVector<Kernel> lbd_;
		DynamicVector<Kernel> ubd_;

		DynamicVector<Kernel> lbx_;
		DynamicVector<Kernel> ubx_;
		DynamicVector<Kernel> lbu_;
		DynamicVector<Kernel> ubu_;

		DynamicMatrix<Kernel> A_;
		DynamicMatrix<Kernel> B_;
		DynamicVector<Kernel> b_;
	};
}
