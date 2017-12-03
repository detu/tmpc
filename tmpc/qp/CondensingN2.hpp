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
	 * \brief CondensingN2 algorithm with O(N^2) runtime in horizon length.
	 * 
	 * Implements the condensing algorithm with quadratic runtime as described in Section 5.4 of Joel's thesis: 
	 * ftp://ftp.esat.kuleuven.be/pub/stadius/ida/reports/13-217.pdf.
	 * 
	 * Manages resources needed for the condensing algorithm of a given size.
	 *
	 * \tparam <Kernel> Kernel for matrix arithmetic.
	 * 
	 * TODO: account for soft constraints.
	 */
	template <typename Kernel_>
	class CondensingN2
	{
	public:
		using Kernel = Kernel_;

		template <typename InIter>
		class Expression
		:	public OcpQpExpressionBase<Expression<InIter>>
		{
			using size_type = typename Kernel::size_t;
			using Matrix = DynamicMatrix<Kernel>;
			using Vector = DynamicVector<Kernel>;
			using Real = typename Kernel::Real;

		public:
			Expression(Expression const&) = delete;
			Expression(Expression &&) = default;

			Expression(CondensingN2& c, InIter first, InIter last)
			:	c_(c)
			,	first_ {first}
			,	last_ {last}
			{
			}


			template <typename QP>
			void implEvalTo(OcpQpBase<QP>& qp) const
			{
				// Initialization of QP variables
				auto nu = first_->size().nu();
				auto nc = first_->size().nc();

				c_.Qc_ = first_->Q();
				submatrix(c_.Rc_, 0, 0, nu, nu) = first_->R();
				//K::left_cols(Sc_, nu) = first_->S();
				submatrix(c_.Sc_, 0, 0, c_.Sc_.rows(), nu) = first_->S();
				c_.qc_ = first_->q();
				//K::head(c_.rc_, nu) = first_->r();
				subvector(c_.rc_, 0, nu) = first_->r();

				//K::top_rows(Cc_, nc) = first_->C();
				submatrix(c_.Cc_, 0, 0, nc, columns(c_.Cc_)) = first_->C();
				//K::top_left_corner(Dc_, nc, nu) = first_->D();
				submatrix(c_.Dc_, 0, 0, nc, nu) = first_->D();
				//K::head(lbd, nc) = first_->lbd();
				subvector(c_.lbd_, 0, nc) = first_->lbd();
				//K::head(ubd, nc) = first_->ubd();
				subvector(c_.ubd_, 0, nc) = first_->ubd();

				//K::head(lbu, nu) = first_->lbu();
				subvector(c_.lbu_, 0, nu) = first_->lbu();
				//K::head(ubu, nu) = first_->ubu();
				subvector(c_.ubu_, 0, nu) = first_->ubu();

				c_.A_ = first_->A();
				c_.B_ = first_->B();
				c_.b_ = first_->b();

				for (auto stage = first_ + 1; stage != last_; ++stage)
				{
					auto const& sz = stage->size();
					// TODO: add nx1() to OcpSize
					auto const nx_next = stage->B().rows();

					// Update Q
					c_.Qc_ += trans(c_.A_) * stage->Q() * c_.A_;

					// Update S
					//K::left_cols(c_.Sc_, nu) += K::trans(A) * stage->Q() * B;
					submatrix(c_.Sc_, 0, 0, c_.Sc_.rows(), nu) += trans(c_.A_) * stage->Q() * c_.B_;
					//K::middle_cols(c_.Sc_, nu, sz.nu()) = K::trans(A) * stage->S();
					submatrix(c_.Sc_, 0, nu, c_.Sc_.rows(), sz.nu()) = trans(c_.A_) * stage->S();

					// Update R
					// This is the hottest line of the algorithm:
					submatrix(c_.Rc_, 0, 0, nu, nu) += trans(c_.B_) * stage->Q() * c_.B_;
					submatrix(c_.Rc_, 0, nu, nu, sz.nu()) = trans(c_.B_) * stage->S();
					submatrix(c_.Rc_, nu, 0, sz.nu(), nu) = trans(stage->S()) * c_.B_;
					submatrix(c_.Rc_, nu, nu, sz.nu(), sz.nu()) = stage->R();

					// Update q
					c_.qc_ += trans(c_.A_) * (stage->Q() * c_.b_ + stage->q());

					// Update r
					//K::head(c_.rc_, nu) += K::trans(B) * (stage->Q() * b + stage->q());
					subvector(c_.rc_, 0, nu) += trans(c_.B_) * (stage->Q() * c_.b_ + stage->q());
					//K::segment(c_.rc_, nu, sz.nu()) = K::trans(stage->S()) * b + stage->r();
					subvector(c_.rc_, nu, sz.nu()) = trans(stage->S()) * c_.b_ + stage->r();

					// Update C
					//K::middle_rows(c_.Cc_, nc, sz.nx() + sz.nc()) << A, stage->C() * A;
					submatrix(c_.Cc_, nc          , 0, sz.nx(), columns(c_.Cc_)) = c_.A_;
					submatrix(c_.Cc_, nc + sz.nx(), 0, sz.nc(), columns(c_.Cc_)) = stage->C() * c_.A_;

					// Update D (which is initialized to 0)
					//K::middle_rows(c_.Dc_, nc, sz.nx() + sz.nc()) << B, K::zero(sz.nx(), cs_.nu() - nu),
					//									  stage->C() * B, stage->D(), K::zero(sz.nc(), cs_.nu() - nu - sz.nu());
					submatrix(c_.Dc_, nc          ,  0, sz.nx(),      nu) = c_.B_;
					submatrix(c_.Dc_, nc + sz.nx(),  0, sz.nc(),      nu) = stage->C() * c_.B_;
					submatrix(c_.Dc_, nc + sz.nx(), nu, sz.nc(), sz.nu()) = stage->D();

					// Update lbd
					//K::segment(lbd, nc, sz.nx() + sz.nc()) << stage->lbx() - b, stage->lbd() - stage->C() * b;
					subvector(c_.lbd_, nc          , sz.nx()) = stage->lbx() - c_.b_;
					subvector(c_.lbd_, nc + sz.nx(), sz.nc()) = stage->lbd() - stage->C() * c_.b_;

					// Update ubd
					//K::segment(ubd, nc, sz.nx() + sz.nc()) << stage->ubx() - b, stage->ubd() - stage->C() * b;
					subvector(c_.ubd_, nc          , sz.nx()) = stage->ubx() - c_.b_;
					subvector(c_.ubd_, nc + sz.nx(), sz.nc()) = stage->ubd() - stage->C() * c_.b_;

					// Uplate lbu, ubu
					subvector(c_.lbu_, nu, sz.nu()) = stage->lbu();
					subvector(c_.ubu_, nu, sz.nu()) = stage->ubu();

					// Update A
					c_.A_ = stage->A() * c_.A_;

					// Update B
					{
						DynamicMatrix<Kernel> B_next(nx_next, nu + sz.nu());
						submatrix(B_next, 0,  0, nx_next,      nu) = stage->A() * c_.B_;
						submatrix(B_next, 0, nu, nx_next, sz.nu()) = stage->B();
						c_.B_ = std::move(B_next);
					}

					// Update b
					c_.b_ = stage->A() * c_.b_ + stage->b();

					// Update indices
					nu += sz.nu();
					nc += sz.nx() + sz.nc();
				}

				// Set the state bounds of the condensed stage
				qp.lbx(first_->lbx());
				qp.ubx(first_->ubx());

				// Assign the output values.				
				qp.A(c_.A_);
				qp.B(c_.B_);
				qp.b(c_.b_);
				qp.C(c_.Cc_);
				qp.D(c_.Dc_);
				qp.lbd(c_.lbd_);
				qp.ubd(c_.ubd_);
				qp.lbu(c_.lbu_);
				qp.ubu(c_.ubu_);
				qp.Q(c_.Qc_);
				qp.R(c_.Rc_);
				qp.S(c_.Sc_);
				qp.q(c_.qc_);
				qp.r(c_.rc_);

				// TODO: implement soft constraints matrices recalculation.
				qp.Zl(sNaN<Real>());
				qp.Zu(sNaN<Real>());
				qp.zl(sNaN<Real>());
				qp.zu(sNaN<Real>());

				// TODO: what happens to the soft constraints index?
			}


			auto const& implSize() const
			{
				return c_.cs_;
			}


		private:
			CondensingN2& c_;
			InIter first_;
			InIter last_;
		};

		template <typename InIter>
		CondensingN2(InIter const& sz_first, InIter const& sz_last)
		:	cs_(condensedOcpSize(sz_first, sz_last))
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
		Expression<InIter> operator()(InIter first, InIter last)
		{
			if (condensedOcpSize(ocpSizeIterator(first), ocpSizeIterator(last)) != cs_)
				throw std::invalid_argument("QP size does not match the condensed size");

			return Expression<InIter> {*this, first, last};
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

		DynamicVector<Kernel> lbu_;
		DynamicVector<Kernel> ubu_;

		DynamicMatrix<Kernel> A_;
		DynamicMatrix<Kernel> B_;
		DynamicVector<Kernel> b_;
	};
}
