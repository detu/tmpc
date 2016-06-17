#pragma once

#include <Eigen/Dense>

namespace tmpc
{
	template<typename MultiStageQP_, typename CondensedQP_>
	void Condense(MultiStageQP_ const& msqp, CondensedQP_& condensed_qp)
	{
		using namespace Eigen;

		if (nIndep(msqp) != condensed_qp.nx())
		{
			std::stringstream msg;
			msg << "Condense(): output CondensedQP has " << condensed_qp.nx()
					<< " variables, but is expected to have " << nIndep(msqp) << ".";
			throw std::invalid_argument(msg.str());
		}

		if (nDep(msqp) + nConstr(msqp) != condensed_qp.nc())
		{
			std::stringstream msg;
			msg << "Condense(): output CondensedQP has " << condensed_qp.nc()
					<< " constraints, but is expected to have " << nDep(msqp) + nConstr(msqp) << ".";
			throw std::invalid_argument(msg.str());
		}

		auto constexpr nX = MultiStageQP_::nX();
		auto constexpr nU = MultiStageQP_::nU();
		auto constexpr nD = MultiStageQP_::nD();
		auto constexpr nDT = MultiStageQP_::nDT();
		auto const nT = msqp.nT();
		auto const n_indep = nIndep(msqp);
		auto constexpr nC = nX + nD;

		MatrixXd M = MatrixXd::Identity(nX, n_indep);
		auto v = Matrix<double, nX, 1>::Zero().eval();

		auto& Hc = condensed_qp.H();
		auto& gc = condensed_qp.g();

		Hc.setZero();
		gc.setZero();

		condensed_qp.lb().template topRows<nX>() = xMin(msqp, 0);
		condensed_qp.ub().template topRows<nX>() = xMax(msqp, 0);

		for (unsigned k = 0; k < nT; ++k)
		{
			auto M_k = M.leftCols(nX + k * nU);

			// Calculate Hessian (H_k) and gradient (g_k) for current time step w.r.t. condensed (independent) variables.
			const auto H_k = msqp.H(k);
			const auto g_k = msqp.g(k);

			const auto Q = H_k.template topLeftCorner<nX, nX>().template selfadjointView<Upper>();
			const auto S = H_k.template topRightCorner<nX, nU>();
			const auto R = H_k.template bottomRightCorner<nU, nU>();

			const auto nn = M_k.cols();
			auto Hc_k = Hc.topLeftCorner(nn + nU, nn + nU);
			Hc_k.topLeftCorner(nn, nn).template triangularView<Upper>() += M_k.transpose() * Q * M_k;
			Hc_k.topRightCorner(nn, nU) += M_k.transpose() * S;
			Hc_k.template bottomRightCorner<nU, nU>().template triangularView<Upper>() += R;

			auto gc_k = gc.topRows(nn + nU);
			gc_k.topRows(nn) += M_k.transpose() * (g_k.template topRows<nX>() + Q * v);
			gc_k.template bottomRows<nU>() += g_k.template bottomRows<nU>() + S.transpose() * v;

			// Set lower and upper bound for independent variables at stage k.
			condensed_qp.lb().template middleRows<nU>(nX + k * nU) = uMin(msqp, k);
			condensed_qp.ub().template middleRows<nU>(nX + k * nU) = uMax(msqp, k);

			// Set path constraints for stage k.
			auto Aconstr_k = condensed_qp.A().template middleRows<nC>(k * nC);
			auto lbAconstr_k = condensed_qp.lbA().template middleRows<nC>(k * nC);
			auto ubAconstr_k = condensed_qp.ubA().template middleRows<nC>(k * nC);
			auto D_k_x = msqp.D(k).template leftCols<nX>();
			auto D_k_u = msqp.D(k).template rightCols<nU>();
			const auto d_ofs = (D_k_x * v).eval();

			Aconstr_k.template topRows<nD>() = D_k_x * M;
			Aconstr_k.template block<nD, nU>(0, nn) = D_k_u;
			lbAconstr_k.template topRows<nD>() = msqp.dMin(k) - d_ofs;
			ubAconstr_k.template topRows<nD>() = msqp.dMax(k) - d_ofs;
			
			// Update M and v.
			const auto A_k = msqp.C(k).template leftCols<nX>();
			const auto B_k = msqp.C(k).template rightCols<nU>();
			M_k = A_k * M_k;
			M.middleCols<nU>(nX + k * nU) = B_k;
			v = A_k * v + msqp.c(k);

			// Set next state bound constraints.
			Aconstr_k.template bottomRows<nX>() = M;

			if (k + 1 < nT)
			{
				lbAconstr_k.template bottomRows<nX>() = xMin(msqp, k + 1) - v;
				ubAconstr_k.template bottomRows<nX>() = xMax(msqp, k + 1) - v;
			}
			else
			{
				lbAconstr_k.template bottomRows<nX>() = msqp.zendMin() - v;
				ubAconstr_k.template bottomRows<nX>() = msqp.zendMax() - v;
			}
		}

		// Set terminal constraints.
		const auto D_k_term = msqp.Dend();
		const auto d_ofs = (D_k_term * v).eval();

		condensed_qp.A().template bottomRows<nDT>() = D_k_term * M;
		condensed_qp.lbA().template bottomRows<nDT>() = msqp.dendMin() - d_ofs;
		condensed_qp.ubA().template bottomRows<nDT>() = msqp.dendMax() - d_ofs;

		// Cost of final state.
		Hc.template triangularView<Upper>() += M.transpose() * msqp.Hend().template selfadjointView<Upper>() * M;
		gc += M.transpose() * (msqp.gend() + msqp.Hend() * v);

		Hc = Hc.template selfadjointView<Upper>();
	}
}
