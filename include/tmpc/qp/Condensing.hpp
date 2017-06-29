#pragma once

#include <tmpc/qp/QpSize.hpp>

#include <blaze/Math.h>

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
	 * \brief Condense a multistage (sparse) QP to a dense QP.
	 *
	 * \tparam <K> Class implementing the Kernel concept
	 * \tparam <D> Class defining problem dimensions
	 * \tparam <MultiStageQP_> Class implementing MultiStageQP concept
	 * \tparam <CondensedQP_> Class implementing DenseQP concept
	 *
	 * TODO: should D be deduced from MultiStageQP_ or not?
	 * */
	template <typename K, typename D, typename MultiStageQP_, typename CondensedQP_>
	void Condense(MultiStageQP_ const& msqp, CondensedQP_& condensed_qp)
	{
		// TODO: condendes_qp should, actually, be a QpStage (?)
		if (condensed_qp.size() != 1 || condensedQpSize(msqp.size()) != condensed_qp.front().size())
		{
			std::stringstream msg;
			msg << "Condense(): output QP has wrong size.";
			throw std::invalid_argument(msg.str());
		}

		auto constexpr nX = D::NX;
		auto constexpr nU = D::NU;
		auto constexpr nD = D::NC;
		auto constexpr nDT = D::NCT;
		auto const nT = msqp.nT();
		auto const n_indep = nIndep(msqp);
		auto constexpr nC = nX + nD;

		typename K::DynamicMatrix M = K::identity(nX, n_indep);
		typename K::template Vector<D::NX> v = K::template zero<D::NX>();

		auto& Hc = condensed_qp.H();
		auto& gc = condensed_qp.g();

		K::set_zero(Hc);
		K::set_zero(gc);

		K::template top_rows<nX>(condensed_qp.lb()) = msqp.get_x_min(0);
		K::template top_rows<nX>(condensed_qp.ub()) = msqp.get_x_max(0);

		for (unsigned k = 0; k < nT; ++k)
		{
			auto M_k = K::left_cols(M, nX + k * nU);

			// Calculate Hessian (H_k) and gradient (g_k) for current time step w.r.t. condensed (independent) variables.
			auto const H_k = get_H(msqp, k);
			auto const g_k = get_g(msqp, k);

			const auto Q = K::selfadjoint_view_upper(K::template top_left_corner    <nX, nX>(H_k));
			const auto S =                           K::template top_right_corner   <nX, nU>(H_k) ;
			const auto R =                           K::template bottom_right_corner<nU, nU>(H_k) ;

			const auto nn = K::cols(M_k);
			auto Hc_k = K::top_left_corner(Hc, nn + nU, nn + nU);
			K::triangular_view_upper(K::top_left_corner(Hc_k, nn, nn)) += K::transpose(M_k) * Q * M_k;
			K::top_right_corner(Hc_k, nn, nU) += K::transpose(M_k) * S;
			K::triangular_view_upper(K::template bottom_right_corner<nU, nU>(Hc_k)) += R;

			auto gc_k = K::top_rows(gc, nn + nU);
			K::top_rows(gc_k, nn) += K::transpose(M_k) * (K::template top_rows<nX>(g_k) + Q * v);
			K::template bottom_rows<nU>(gc_k) += K::template bottom_rows<nU>(g_k) + K::transpose(S) * v;

			// Set lower and upper bound for independent variables at stage k.
			K::template middle_rows<nU>(condensed_qp.lb(), nX + k * nU) = msqp.get_u_min(k);
			K::template middle_rows<nU>(condensed_qp.ub(), nX + k * nU) = msqp.get_u_max(k);

			// Set path constraints for stage k.
			auto Aconstr_k   = K::template middle_rows<nC>(condensed_qp.A  (), k * nC);
			auto lbAconstr_k = K::template middle_rows<nC>(condensed_qp.lbA(), k * nC);
			auto ubAconstr_k = K::template middle_rows<nC>(condensed_qp.ubA(), k * nC);
			auto D_k_x = msqp.get_C(k);
			auto D_k_u = msqp.get_D(k);
			const auto d_ofs = K::eval(D_k_x * v);

			K::template top_rows<nD>(Aconstr_k) = D_k_x * M;
			K::template block<nD, nU>(Aconstr_k, 0, nn) = D_k_u;
			K::template top_rows<nD>(lbAconstr_k) = msqp.get_d_min(k) - d_ofs;
			K::template top_rows<nD>(ubAconstr_k) = msqp.get_d_max(k) - d_ofs;
			
			// Update M and v.
			const auto A_k = msqp.get_A(k);
			const auto B_k = msqp.get_B(k);
			M_k = A_k * M_k;
			K::template middle_cols<nU>(M, nX + k * nU) = B_k;
			v = A_k * v + msqp.get_b(k);

			// Set next state bound constraints.
			K::template bottom_rows<nX>(Aconstr_k) = M;

			if (k + 1 < nT)
			{
				K::template bottom_rows<nX>(lbAconstr_k) = msqp.get_x_min(k + 1) - v;
				K::template bottom_rows<nX>(ubAconstr_k) = msqp.get_x_max(k + 1) - v;
			}
			else
			{
				K::template bottom_rows<nX>(lbAconstr_k) = msqp.get_x_min(nT) - v;
				K::template bottom_rows<nX>(ubAconstr_k) = msqp.get_x_max(nT) - v;
			}
		}

		// Set terminal constraints.
		const auto D_k_term = msqp.get_C_end();
		const auto d_ofs = K::eval(D_k_term * v);

		K::template bottom_rows<nDT>(condensed_qp.A  ()) = D_k_term * M;
		K::template bottom_rows<nDT>(condensed_qp.lbA()) = msqp.get_d_end_min() - d_ofs;
		K::template bottom_rows<nDT>(condensed_qp.ubA()) = msqp.get_d_end_max() - d_ofs;

		// Cost of final state.
		K::triangular_view_upper(Hc) += K::transpose(M) * K::selfadjoint_view_upper(msqp.get_Q(nT)) * M;
		gc += K::transpose(M) * (msqp.get_q(nT) + msqp.get_Q(nT) * v);

		Hc = K::selfadjoint_view_upper(Hc);
	}

	/**
	 * \brief Condense a multistage QP to a dense QP.
	 *
	 * \tparam <K> Class implementing the Kernel concept
	 * \tparam <InIter> Class defining input iterator to multistage QP stages
	 * \tparam <CondensedQP_> Class implementing a QpStage concept
	 *
	 * TODO: should D be deduced from MultiStageQP_ or not?
	 * */
	template <typename K, typename InIter, typename CondensedQP_>
	void Condense(InIter qp_begin, InIter qp_end, CondensedQP_& condensed_qp)
	{
		if (qp_begin == qp_end)
			throw std::invalid_argument("Condense(): the input QpStage must not be empty");

		QpSize const cs = condensedQpSize(tmpc::qpSizeIterator(qp_begin), tmpc::qpSizeIterator(qp_end));

		typename K::DynamicMatrix Qc(cs.nx(), cs.nx());
		typename K::DynamicMatrix Rc(cs.nu(), cs.nu());
		typename K::DynamicMatrix Sc(cs.nx(), cs.nu());
		typename K::DynamicVector qc(cs.nx());
		typename K::DynamicVector rc(cs.nu());

		typename K::DynamicMatrix Cc(cs.nc(), cs.nx());
		typename K::DynamicMatrix Dc(cs.nc(), cs.nu());
		typename K::DynamicVector lbd(cs.nc());
		typename K::DynamicVector ubd(cs.nc());

		typename K::DynamicVector lbu(cs.nu());
		typename K::DynamicVector ubu(cs.nu());

		// Initialization of QP variables
		auto nu = qp_begin->size().nu();
		auto nc = qp_begin->size().nc();

		Qc = qp_begin->get_Q();
		K::top_left_corner(Rc, nu, nu) = qp_begin->get_R();
		K::left_cols(Sc, nu) = qp_begin->get_S();
		qc = qp_begin->get_q();
		K::head(rc, nu) = qp_begin->get_r();

		K::top_rows(Cc, nc) = qp_begin->get_C();
		K::top_left_corner(Dc, nc, nu) = qp_begin->get_D();
		K::head(lbd, nc) = qp_begin->get_lbd();
		K::head(ubd, nc) = qp_begin->get_ubd();

		K::head(lbu, nu) = qp_begin->get_lbu();
		K::head(ubu, nu) = qp_begin->get_ubu();

		typename K::DynamicMatrix A = qp_begin->get_A();
		typename K::DynamicMatrix B = qp_begin->get_B();
		typename K::DynamicVector b = qp_begin->get_b();

		for (auto stage = qp_begin + 1; stage != qp_end; ++stage)
		{
			auto const& sz = stage->size();
			// TODO: add nx1() to QpSize
			auto const nx_next = K::rows(stage->get_B());

			// Update Q
			Qc += K::transpose(A) * stage->get_Q() * A;

			// Update S
			K::left_cols(Sc, nu) += K::transpose(A) * stage->get_Q() * B;
			K::middle_cols(Sc, nu, sz.nu()) = K::transpose(A) * stage->get_S();

			// Update R
			K::top_left_corner(Rc, nu, nu) += K::transpose(B) * stage->get_Q() * B;
			K::block(Rc, 0, nu, nu, sz.nu()) = K::transpose(B) * stage->get_S();
			K::block(Rc, nu, 0, sz.nu(), nu) = K::transpose(stage->get_S()) * B;
			K::block(Rc, nu, nu, sz.nu(), sz.nu()) = stage->get_R();

			// Update q
			qc += K::transpose(A) * (stage->get_Q() * b + stage->get_q());

			// Update r
			K::head(rc, nu) += K::transpose(B) * (stage->get_Q() * b + stage->get_q());
			K::segment(rc, nu, sz.nu()) = K::transpose(stage->get_S()) * b + stage->get_r();

			// Update C
			K::middle_rows(Cc, nc, sz.nx() + sz.nc()) << A, stage->get_C() * A;

			// Update D
			K::middle_rows(Dc, nc, sz.nx() + sz.nc()) << B, K::zero(sz.nx(), cs.nu() - nu),
			                                      stage->get_C() * B, stage->get_D(), K::zero(sz.nc(), cs.nu() - nu - sz.nu());

			// Update lbd
			K::segment(lbd, nc, sz.nx() + sz.nc()) << stage->get_lbx() - b, stage->get_lbd() - stage->get_C() * b;

			// Update ubd
			K::segment(ubd, nc, sz.nx() + sz.nc()) << stage->get_ubx() - b, stage->get_ubd() - stage->get_C() * b;

			// Update A
			{
				typename K::DynamicMatrix A_next = stage->get_A() * A;
				A.resizeLike(A_next);
				A.swap(A_next);
			}

			// Update B
			{
				typename K::DynamicMatrix B_next(nx_next, nu + sz.nu());
				B_next << stage->get_A() * B, stage->get_B();
				B.resizeLike(B_next);
				B.swap(B_next);
			}

			// Update b
			{
				typename K::DynamicMatrix b_next = stage->get_A() * b + stage->get_b();
				b.resizeLike(b_next);
				b.swap(b_next);
			}

			// Update indices
			nu += sz.nu();
			nc += sz.nx() + sz.nc();
		}

		// Fill the condensed stage
		condensed_qp.set_Q(Qc);
		condensed_qp.set_R(Rc);
		condensed_qp.set_S(Sc);
		condensed_qp.set_q(qc);
		condensed_qp.set_r(rc);

		condensed_qp.set_C(Cc);
		condensed_qp.set_D(Dc);
		condensed_qp.set_lbd(lbd);
		condensed_qp.set_ubd(ubd);

		condensed_qp.set_A(A);
		condensed_qp.set_B(B);
		condensed_qp.set_b(b);

		condensed_qp.set_lbx(qp_begin->get_lbx());
		condensed_qp.set_ubx(qp_begin->get_ubx());
		condensed_qp.set_lbu(lbu);
		condensed_qp.set_ubu(ubu);
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

		blaze::DynamicMatrix<Scalar> Qc(cs.nx(), cs.nx());
		blaze::DynamicMatrix<Scalar> Rc(cs.nu(), cs.nu());
		blaze::DynamicMatrix<Scalar> Sc(cs.nx(), cs.nu());
		blaze::DynamicVector<Scalar> qc(cs.nx());
		blaze::DynamicVector<Scalar> rc(cs.nu());

		blaze::DynamicMatrix<Scalar> Cc(cs.nc(), cs.nx());
		blaze::DynamicMatrix<Scalar> Dc(cs.nc(), cs.nu(), 0);
		blaze::DynamicVector<Scalar> lbd(cs.nc());
		blaze::DynamicVector<Scalar> ubd(cs.nc());

		blaze::DynamicVector<Scalar> lbu(cs.nu());
		blaze::DynamicVector<Scalar> ubu(cs.nu());

		// Initialization of QP variables
		auto nu = first->size().nu();
		auto nc = first->size().nc();

		Qc = first->Q();
		submatrix(Rc, 0, 0, nu, nu) = first->R();
		submatrix(Sc, 0, 0, rows(Sc), nu) = first->S();
		qc = first->q();
		subvector(rc, 0, nu) = first->r();

		submatrix(Cc, 0, 0, nc, columns(Cc)) = first->C();
		submatrix(Dc, 0, 0, nc, nu) = first->D();
		subvector(lbd, 0, nc) = first->lbd();
		subvector(ubd, 0, nc) = first->ubd();

		subvector(lbu, 0, nu) = first->lbu();
		subvector(ubu, 0, nu) = first->ubu();

		blaze::DynamicMatrix<Scalar> A = first->A();
		blaze::DynamicMatrix<Scalar> B = first->B();
		blaze::DynamicVector<Scalar> b = first->b();

		for (auto stage = first + 1; stage != last; ++stage)
		{
			auto const& sz = stage->size();
			// TODO: add nx1() to QpSize
			auto const nx_next = rows(stage->B());

			// Update Q
			Qc += trans(A) * stage->Q() * A;

			// Update S
			submatrix(Sc, 0, 0, rows(Sc), nu) += trans(A) * stage->Q() * B;
			submatrix(Sc, 0, nu, rows(Sc), sz.nu()) = trans(A) * stage->S();

			// Update R
			submatrix(Rc, 0, 0, nu, nu) += trans(B) * stage->Q() * B;
			submatrix(Rc, 0, nu, nu, sz.nu()) = trans(B) * stage->S();
			submatrix(Rc, nu, 0, sz.nu(), nu) = trans(stage->S()) * B;
			submatrix(Rc, nu, nu, sz.nu(), sz.nu()) = stage->R();

			// Update q
			qc += trans(A) * (stage->Q() * b + stage->q());

			// Update r
			subvector(rc, 0, nu) += trans(B) * (stage->Q() * b + stage->q());
			subvector(rc, nu, sz.nu()) = trans(stage->S()) * b + stage->r();

			// Update C
			submatrix(Cc, nc          , 0, sz.nx(), columns(Cc)) = A;
			submatrix(Cc, nc + sz.nx(), 0, sz.nc(), columns(Cc)) = stage->C() * A;

			// Update D (which is initialized to 0)
			submatrix(Dc, nc          ,  0, sz.nx(),      nu) = B;
			submatrix(Dc, nc + sz.nx(),  0, sz.nc(),      nu) = stage->C() * B;
			submatrix(Dc, nc + sz.nx(), nu, sz.nc(), sz.nu()) = stage->D();

			// Update lbd
			subvector(lbd, nc          , sz.nx()) = stage->lbx() - b;
			subvector(lbd, nc + sz.nx(), sz.nc()) = stage->lbd() - stage->C() * b;

			// Update ubd
			subvector(ubd, nc          , sz.nx()) = stage->ubx() - b;
			subvector(ubd, nc + sz.nx(), sz.nc()) = stage->ubd() - stage->C() * b;

			// Update A
			A = stage->A() * A;

			// Update B
			{
				blaze::DynamicMatrix<Scalar> B_next(nx_next, nu + sz.nu());
				submatrix(B_next, 0,  0, nx_next,      nu) = stage->A() * B;
				submatrix(B_next, 0, nu, nx_next, sz.nu()) = stage->B();
				B = std::move(B_next);
			}

			// Update b
			b = stage->A() * b + stage->b();

			// Update indices
			nu += sz.nu();
			nc += sz.nx() + sz.nc();
		}

		// Fill the condensed stage
		condensed.Q() = Qc;
		condensed.R() = Rc;
		condensed.S() = Sc;
		condensed.q() = qc;
		condensed.r() = rc;

		condensed.C() = Cc;
		condensed.D() = Dc;
		condensed.lbd() = lbd;
		condensed.ubd() = ubd;

		condensed.A() = A;
		condensed.B() = B;
		condensed.b() = b;

		condensed.lbx() = first->lbx();
		condensed.ubx() = first->ubx();
		condensed.lbu() = lbu;
		condensed.ubu() = ubu;
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
	CondenseExpression<Scalar, InIter> condenseBlaze(InIter first, InIter last)
	{
		return CondenseExpression<Scalar, InIter>(first, last);
	}
}
