#pragma once

#include <tmpc/qp/QpSize.hpp>
#include <tmpc/Matrix.hpp>

#include <vector>
#include <sstream>
#include <stdexcept>
#include <numeric>

namespace tmpc
{
	/**
	 * Resulting QP size after full condensing.
	 */
	QpSize condensedQpSize(std::vector<QpSize> const& sz);

	template <typename InIter>
	QpSize condensedQpSize(InIter sz_begin, InIter sz_end)
	{
		if (sz_begin == sz_end)
			throw std::invalid_argument("condensedQpSize(): QpSize range must be not empty");

		return QpSize(
			sz_begin->nx(),
			std::accumulate(sz_begin, sz_end, std::size_t{0},
				[] (std::size_t n, QpSize const& s) { return n + s.nu(); }),
			std::accumulate(sz_begin, sz_end, std::size_t{0},
				[] (std::size_t n, QpSize const& s) { return n + s.nc(); })
			+ std::accumulate(sz_begin + 1, sz_end, std::size_t{0},
					[] (std::size_t n, QpSize const& s) { return n + s.nx(); })
			);
	}

	/**
	 * \brief Condensing algorithm.
	 *
	 * Manages resources needed for the condensing algorithm of a given size.
	 *
	 * \tparam <Scalar> Scalar type used for intermediate variables in the condensing algorithm.
	 */
	template <typename Scalar>
	class Condensing
	{
	public:
		template <typename InIter>
		Condensing(InIter const& sz_first, InIter const& sz_last);
	};

	/**
	 * \brief An expression defining a call to condensing algorithm.
	 * \tparam <Scalar> Scalar type used for intermediate variables in the condensing algorithm.
	 * \tparam <InIter> Class defining input iterator to multistage QP stages
	 */
	template <typename Scalar, typename InIter>
	class CondenseExpression
	{
	public:
		CondenseExpression(InIter const& first, InIter const& last)
		:	first_(first)
		,	last_(last)
		{
			if (first == last)
				throw std::invalid_argument("CondenseExpression(): the input range must not be empty");
		}

		InIter const& first() const { return first_; }
		InIter const& last() const { return last_; }

	private:
		InIter const first_;
		InIter const last_;
	};

	/**
	 * \brief Performs the actual condensing.
	 *
	 * \tparam <CondensedStage> Class implementing a QpStage concept
	 */
	template <typename CondensedStage, typename Scalar, typename InIter>
	void assign(CondensedStage& condensed, CondenseExpression<Scalar, InIter> const& rhs)
	{
		auto const& first = rhs.first();
		auto const& last = rhs.last();
		QpSize const cs = condensedQpSize(tmpc::qpSizeIterator(first), tmpc::qpSizeIterator(last));

		DynamicMatrix<Scalar> Qc(cs.nx(), cs.nx());
		DynamicMatrix<Scalar> Rc(cs.nu(), cs.nu());
		DynamicMatrix<Scalar> Sc(cs.nx(), cs.nu());
		DynamicVector<Scalar> qc(cs.nx());
		DynamicVector<Scalar> rc(cs.nu());

		DynamicMatrix<Scalar> Cc(cs.nc(), cs.nx());
		DynamicMatrix<Scalar> Dc(cs.nc(), cs.nu(), 0);
		DynamicVector<Scalar> lbd(cs.nc());
		DynamicVector<Scalar> ubd(cs.nc());

		DynamicVector<Scalar> lbu(cs.nu());
		DynamicVector<Scalar> ubu(cs.nu());

		// Initialization of QP variables
		auto nu = first->size().nu();
		auto nc = first->size().nc();

		Qc = first->get_Q();
		submatrix(Rc, 0, 0, nu, nu) = first->get_R();
		//K::left_cols(Sc, nu) = first_->get_S();
		submatrix(Sc, 0, 0, Sc.rows(), nu) = first->get_S();
		qc = first->get_q();
		//K::head(rc, nu) = first_->get_r();
		subvector(rc, 0, nu) = first->get_r();

		//K::top_rows(Cc, nc) = first_->get_C();
		submatrix(Cc, 0, 0, nc, columns(Cc)) = first->get_C();
		//K::top_left_corner(Dc, nc, nu) = first_->get_D();
		submatrix(Dc, 0, 0, nc, nu) = first->get_D();
		//K::head(lbd, nc) = first_->get_lbd();
		subvector(lbd, 0, nc) = first->get_lbd();
		//K::head(ubd, nc) = first_->get_ubd();
		subvector(ubd, 0, nc) = first->get_ubd();

		//K::head(lbu, nu) = first_->get_lbu();
		subvector(lbu, 0, nu) = first->get_lbu();
		//K::head(ubu, nu) = first_->get_ubu();
		subvector(ubu, 0, nu) = first->get_ubu();

		DynamicMatrix<Scalar> A = first->get_A();
		DynamicMatrix<Scalar> B = first->get_B();
		DynamicVector<Scalar> b = first->get_b();

		for (auto stage = first + 1; stage != last; ++stage)
		{
			auto const& sz = stage->size();
			// TODO: add nx1() to QpSize
			auto const nx_next = stage->get_B().rows();

			// Update Q
			Qc += trans(A) * stage->get_Q() * A;

			// Update S
			//K::left_cols(Sc, nu) += K::trans(A) * stage->get_Q() * B;
			submatrix(Sc, 0, 0, Sc.rows(), nu) += trans(A) * stage->get_Q() * B;
			//K::middle_cols(Sc, nu, sz.nu()) = K::trans(A) * stage->get_S();
			submatrix(Sc, 0, nu, Sc.rows(), sz.nu()) = trans(A) * stage->get_S();

			// Update R
			//K::top_left_corner(Rc, nu, nu) += K::trans(B) * stage->get_Q() * B;
			submatrix(Rc, 0, 0, nu, nu) += trans(B) * stage->get_Q() * B;
			//K::block(Rc, 0, nu, nu, sz.nu()) = K::trans(B) * stage->get_S();
			submatrix(Rc, 0, nu, nu, sz.nu()) = trans(B) * stage->get_S();
			//K::block(Rc, nu, 0, sz.nu(), nu) = K::trans(stage->get_S()) * B;
			submatrix(Rc, nu, 0, sz.nu(), nu) = trans(stage->get_S()) * B;
			//K::block(Rc, nu, nu, sz.nu(), sz.nu()) = stage->get_R();
			submatrix(Rc, nu, nu, sz.nu(), sz.nu()) = stage->get_R();

			// Update q
			qc += trans(A) * (stage->get_Q() * b + stage->get_q());

			// Update r
			//K::head(rc, nu) += K::trans(B) * (stage->get_Q() * b + stage->get_q());
			subvector(rc, 0, nu) += trans(B) * (stage->get_Q() * b + stage->get_q());
			//K::segment(rc, nu, sz.nu()) = K::trans(stage->get_S()) * b + stage->get_r();
			subvector(rc, nu, sz.nu()) = trans(stage->get_S()) * b + stage->get_r();

			// Update C
			//K::middle_rows(Cc, nc, sz.nx() + sz.nc()) << A, stage->get_C() * A;
			submatrix(Cc, nc          , 0, sz.nx(), columns(Cc)) = A;
			submatrix(Cc, nc + sz.nx(), 0, sz.nc(), columns(Cc)) = stage->get_C() * A;

			// Update D (which is initialized to 0)
			//K::middle_rows(Dc, nc, sz.nx() + sz.nc()) << B, K::zero(sz.nx(), cs.nu() - nu),
			//									  stage->get_C() * B, stage->get_D(), K::zero(sz.nc(), cs.nu() - nu - sz.nu());
			submatrix(Dc, nc          ,  0, sz.nx(),      nu) = B;
			submatrix(Dc, nc + sz.nx(),  0, sz.nc(),      nu) = stage->get_C() * B;
			submatrix(Dc, nc + sz.nx(), nu, sz.nc(), sz.nu()) = stage->get_D();

			// Update lbd
			//K::segment(lbd, nc, sz.nx() + sz.nc()) << stage->get_lbx() - b, stage->get_lbd() - stage->get_C() * b;
			subvector(lbd, nc          , sz.nx()) = stage->get_lbx() - b;
			subvector(lbd, nc + sz.nx(), sz.nc()) = stage->get_lbd() - stage->get_C() * b;

			// Update ubd
			//K::segment(ubd, nc, sz.nx() + sz.nc()) << stage->get_ubx() - b, stage->get_ubd() - stage->get_C() * b;
			subvector(ubd, nc          , sz.nx()) = stage->get_ubx() - b;
			subvector(ubd, nc + sz.nx(), sz.nc()) = stage->get_ubd() - stage->get_C() * b;

			// Update A
			A = stage->get_A() * A;

			// Update B
			{
				DynamicMatrix<Scalar> B_next(nx_next, nu + sz.nu());
				submatrix(B_next, 0,  0, nx_next,      nu) = stage->get_A() * B;
				submatrix(B_next, 0, nu, nx_next, sz.nu()) = stage->get_B();
				B = std::move(B_next);
			}

			// Update b
			b = stage->get_A() * b + stage->get_b();

			// Update indices
			nu += sz.nu();
			nc += sz.nx() + sz.nc();
		}

		// Fill the condensed stage
		condensed.set_Q(Qc);
		condensed.set_R(Rc);
		condensed.set_S(Sc);
		condensed.set_q(qc);
		condensed.set_r(rc);

		condensed.set_C(Cc);
		condensed.set_D(Dc);
		condensed.set_lbd(lbd);
		condensed.set_ubd(ubd);

		condensed.set_A(A);
		condensed.set_B(B);
		condensed.set_b(b);

		condensed.set_lbx(first->get_lbx());
		condensed.set_ubx(first->get_ubx());
		condensed.set_lbu(lbu);
		condensed.set_ubu(ubu);
	}

	/**
	 * \brief Lazy expression for condensing of a multistage QP to a dense QP.
	 *
	 * This is a version using Blaze.
	 *
	 * \tparam <Scalar> Scalar type used for intermediate variables in the condensing algorithm.
	 * \tparam <InIter> Class defining input iterator to multistage QP stages
	 *
	 * TODO: should D be deduced from MultiStageQP_ or not?
	 * */
	template <typename Scalar, typename InIter>
	CondenseExpression<Scalar, InIter> condense(InIter first, InIter last)
	{
		return CondenseExpression<Scalar, InIter>(first, last);
	}
}
