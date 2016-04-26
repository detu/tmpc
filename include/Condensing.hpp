#pragma once

#include <Eigen/Dense>

namespace camels
{
	template<typename MultiStageQP_, typename CondensedQP_>
	void Condense(MultiStageQP_ const& msqp, CondensedQP_& condensed_qp)
	{
		using namespace Eigen;

		if (msqp.nIndep() != condensed_qp.nx())
		{
			std::stringstream msg;
			msg << "Condense(): output CondensedQP has " << condensed_qp.nx()
					<< " variables, but is expected to have " << msqp.nIndep() << ".";
			throw std::invalid_argument(msg.str());
		}

		if (msqp.nDep() + msqp.nConstr() != condensed_qp.nc())
		{
			std::stringstream msg;
			msg << "Condense(): output CondensedQP has " << condensed_qp.nc()
					<< " constraints, but is expected to have " << msqp.nDep() + msqp.nConstr() << ".";
			throw std::invalid_argument(msg.str());
		}

		auto const nX = msqp.nX();
		auto const nU = msqp.nU();
		auto const nD = msqp.nD();
		auto const nDT = msqp.nDT();
		auto const nT = msqp.nT();
		auto const nIndep = msqp.nIndep();
		auto const nC = nX + nD;

		MatrixXd M = MatrixXd::Identity(nX, nIndep);
		VectorXd v = VectorXd::Zero(nX);

		auto& Hc = condensed_qp.H();
		auto& gc = condensed_qp.g();

		Hc.setZero();
		gc.setZero();

		condensed_qp.lb().topRows(nX) = msqp.xMin(0);
		condensed_qp.ub().topRows(nX) = msqp.xMax(0);

		for (unsigned k = 0; k < nT; ++k)
		{
			auto M_k = M.leftCols(nX + k * nU);

			// Calculate Hessian (H_k) and gradient (g_k) for current time step w.r.t. condensed (independent) variables.
			const auto H_k = msqp.H(k);
			const auto g_k = msqp.g(k);

			const auto Q = H_k.topLeftCorner(nX, nX).template selfadjointView<Upper>();
			const auto S = H_k.topRightCorner(nX, nU);
			const auto R = H_k.bottomRightCorner(nU, nU);

			const auto nn = M_k.cols();
			auto Hc_k = Hc.topLeftCorner(nn + nU, nn + nU);
			Hc_k.topLeftCorner(nn, nn).template triangularView<Upper>() += M_k.transpose() * Q * M_k;
			Hc_k.topRightCorner(nn, nU) += M_k.transpose() * S;
			Hc_k.bottomRightCorner(nU, nU).template triangularView<Upper>() += R;

			auto gc_k = gc.topRows(nn + nU);
			gc_k.topRows(nn) += M_k.transpose() * (g_k.topRows(nX) + Q * v);
			gc_k.bottomRows(nU) += g_k.bottomRows(nU) + S.transpose() * v;

			// Set lower and upper bound for independent variables at stage k.
			condensed_qp.lb().middleRows(nX + k * nU, nU) = msqp.uMin(k);
			condensed_qp.ub().middleRows(nX + k * nU, nU) = msqp.uMax(k);

			// Set path constraints for stage k.
			auto Aconstr_k = condensed_qp.A().middleRows(k * nC, nC);
			auto lbAconstr_k = condensed_qp.lbA().middleRows(k * nC, nC);
			auto ubAconstr_k = condensed_qp.ubA().middleRows(k * nC, nC);
			auto D_k_x = msqp.D(k).leftCols(nX);
			auto D_k_u = msqp.D(k).rightCols(nU);
			const auto d_ofs = (D_k_x * v).eval();

			Aconstr_k.topRows(nD) = D_k_x * M;
			Aconstr_k.block(0, nn, nD, nU) = D_k_u;
			lbAconstr_k.topRows(nD) = msqp.dMin(k) - d_ofs;
			ubAconstr_k.topRows(nD) = msqp.dMax(k) - d_ofs;
			
			// Update M and v.
			const auto A_k = msqp.C(k).leftCols(nX);
			const auto B_k = msqp.C(k).rightCols(nU);
			M_k = A_k * M_k;
			M.middleCols(nX + k * nU, nU) = B_k;
			v = A_k * v + msqp.c(k);

			// Set next state bound constraints.
			Aconstr_k.bottomRows(nX) = M;
			lbAconstr_k.bottomRows(nX) = msqp.xMin(k + 1) - v;
			ubAconstr_k.bottomRows(nX) = msqp.xMax(k + 1) - v;
		}

		// Set terminal constraints.
		const auto D_k_term = msqp.D(nT);
		const auto d_ofs = (D_k_term * v).eval();

		condensed_qp.A().bottomRows(nDT) = D_k_term * M;
		condensed_qp.lbA().bottomRows(nDT) = msqp.dMin(nT) - d_ofs;
		condensed_qp.ubA().bottomRows(nDT) = msqp.dMax(nT) - d_ofs;

		// Cost of final state.
		Hc.template triangularView<Upper>() += M.transpose() * msqp.H(nT).template selfadjointView<Upper>() * M;
		gc += M.transpose() * (msqp.g(nT) + msqp.H(nT) * v);

		Hc = Hc.template selfadjointView<Upper>();
	}
}
