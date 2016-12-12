#pragma once

#include <tmpc/qp/QpSize.hpp>

#include <vector>
#include <sstream>
#include <stdexcept>

namespace tmpc
{
	/**
	 * Resulting QP size after full condensing.
	 */
	QpSize CondensedQpSize(std::vector<QpSize> const& sz);

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
		if (condensed_qp.size().size() != 1 || CondensedQpSize(msqp.size()) != condensed_qp.size().front())
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
}
