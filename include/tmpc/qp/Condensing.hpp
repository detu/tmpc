#pragma once

#include <tmpc/qp/QpSize.hpp>

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
	template <typename K, typename InIter, typename MultiStageQP_, typename CondensedQP_>
	void Condense(InIter qp_begin, InIter qp_end, CondensedQP_& condensed_qp)
	{
		if (qp_begin == qp_end)
			throw std::invalid_argument("Condense(): the input QpStage must be not empty");

		QpSize const condensed_size = CondensedQpSize(tmpc::qpSizeIterator(qp_begin), tmpc::qpSizeIterator(qp_end));
		if (condensed_size != condensed_qp.size())
		{
			std::stringstream msg;
			msg << "Condense(): output QP stage has wrong size.";
			throw std::invalid_argument(msg.str());
		}

		//auto constexpr nX = D::NX;
		//auto constexpr nU = D::NU;
		//auto constexpr nD = D::NC;
		//auto constexpr nDT = D::NCT;
		//auto constexpr nC = nX + nD;

		auto const nX = condensed_size.nx();
		auto const N = std::distance(qp_begin, qp_end) - 1;
		auto const n_indep = condensed_size.nx() + condensed_size.nu();

		K::DynamicMatrix M_k = K::identity(condensed_size.nx(), condensed_size.nx());
		K::DynamicVector v = K::zero(nX);

		K::DynamicMatrix Hc(n_indep, n_indep);
		K::DynamicVector gc(n_indep);
		//K::DynamicMatrix Qc(condensed_size.nx(), condensed_size.nx());
		K::DynamicMatrix A = K::zero(condensed_size.nc(), condensed_size.nx());
		K::DynamicVector lbA(condensed_size.nc()), ubA(condensed_size.nc());

		K::set_zero(Hc);
		K::set_zero(gc);

		K::DynamicVector lbu(condensed_size.nu()), ubu(condensed_size.nu());

		condensed_qp.set_lbx(qp_begin->get_x_min());
		condensed_qp.set_ubx(qp_begin->get_x_max());

		auto stage = qp_begin;
		std::size_t nn = condensed_size.nx();
		std::size_t constr_idx = 0;

		for (unsigned k = 0; k < N; ++k, ++stage)
		{
			auto const& sz = stage->size();
			//auto M_k = K::left_cols(M, nX + k * nU);

			// Calculate Hessian (H_k) and gradient (g_k) for current time step w.r.t. condensed (independent) variables.
			//auto const H_k = get_H(msqp, k);
			//auto const g_k = get_g(msqp, k);

			const auto Q = K::selfadjoint_view_upper(stage->get_Q());
			const auto S =                           stage->get_S() ;
			const auto R =                           stage->get_R() ;

			auto const nU = sz.nu();
			auto const nD = sz.nc();
			auto const nC = sz.nx() + sz.nc();

			auto Hc_k = K::top_left_corner(Hc, nn + nU, nn + nU);
			K::triangular_view_upper(K::top_left_corner(Hc_k, nn, nn)) += K::transpose(M_k) * Q * M_k;
			K::top_right_corner(Hc_k, nn, nU) += K::transpose(M_k) * S;
			K::triangular_view_upper(K::bottom_right_corner(Hc_k, nU, nU)) += R;

			auto gc_k = K::top_rows(gc, nn + nU);
			K::top_rows(gc_k, nn) += K::transpose(M_k) * (stage->get_q() + Q * v);
			K::bottom_rows(gc_k, nU) += stage->get_r() + K::transpose(S) * v;

			// Set lower and upper bound for independent variables at stage k.
			K::middle_rows(lbu, nn - nX, nU) = stage->get_u_min();
			K::middle_rows(ubu, nn - nX, nU) = stage->get_u_max();

			// Set path constraints for stage k.
			auto Aconstr_k   = K::middle_rows(A  , constr_idx, nC);
			auto lbAconstr_k = K::middle_rows(lbA, constr_idx, nC);
			auto ubAconstr_k = K::middle_rows(ubA, constr_idx, nC);
			auto D_k_x = stage->get_C();
			auto D_k_u = stage->get_D();
			const auto d_ofs = K::eval(D_k_x * v);

			K::top_left_corner(Aconstr_k, nD, nn) = D_k_x * M_k;
			K::block(Aconstr_k, 0, nn, nD, nU) = D_k_u;
			K::top_rows(lbAconstr_k, nD) = stage->get_lbd() - d_ofs;
			K::top_rows(ubAconstr_k, nD) = stage->get_ubd() - d_ofs;

			// Update M and v.
			const auto A_k = stage->get_A();
			const auto B_k = stage->get_B();

			{
				K::DynamicMatrix M_next(A_k.rows(), M_k.cols() + B_k.cols());
				M_next << A_k * M_k, B_k;
				M_k.swap(M_next);
			}

			v = A_k * v + stage->get_b();

			if (k < nT)
			{
				// Set next state bound constraints.
				K::bottom_left_corner(Aconstr_k, sz.nx(), nn) = M_k;
				K::template bottom_rows<nX>(lbAconstr_k) = msqp.get_x_min(k + 1) - v;
				K::template bottom_rows<nX>(ubAconstr_k) = msqp.get_x_max(k + 1) - v;
			}

			nn += nU;
			constr_idx += nC;
		}

		// Set terminal constraints.
		const auto D_k_term = msqp.get_C_end();
		const auto d_ofs = K::eval(D_k_term * v);

		K::template bottom_rows<nDT>(condensed_qp.A  ()) = D_k_term * M;
		K::template bottom_rows<nDT>(condensed_qp.lbA()) = msqp.get_d_end_min() - d_ofs;
		K::template bottom_rows<nDT>(condensed_qp.ubA()) = msqp.get_d_end_max() - d_ofs;

		// Cost of final state.
		K::triangular_view_upper(Hc) += K::transpose(M) * K::selfadjoint_view_upper(msqp.get_Q(nT)) * M;

		condensed_qp.set_g(gc + K::transpose(M) * (msqp.get_q(nT) + msqp.get_Q(nT) * v));
		condensed_qp.set_H(K::selfadjoint_view_upper(Hc));
		condensed_qp.set_lbu(lbu);
		condensed_qp.set_ubu(ubu);
	}
}
