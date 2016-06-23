#pragma once

#include <stdexcept>

namespace tmpc_test
{
	namespace qp_problems
	{
		template <typename Problem>
		void problem_0(Problem& qp)
		{
			unsigned const NX = 2;
			unsigned const NU = 1;
			unsigned const NC = 0;
			unsigned const NCT = 0;
			unsigned const NT = 2;

			static_assert(qp.nX() == NX && qp.nU() == NU && qp.nD() == NC && qp.nDT() == NCT, "Problem sizes must match");

			if (qp.nT() != NT)
				throw std::logic_error("Invalid size of the QP problem");

			set_zMin(qp, 0, -1.);	set_zMax(qp, 0, 1.);
			set_zMin(qp, 1, -1.);	set_zMax(qp, 1, 1.);
			set_zendMin(qp, -1.);	set_zendMax(qp, 1.);

			// Stage 0
			typename Problem::StageHessianMatrix H0;
			H0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
			H0 = H0.transpose() * H0;	// Make positive definite.

			typename Problem::StageGradientVector g0;
			g0 << 0., 0., 0.;

			const Eigen::MatrixXd Q0 = H0.topLeftCorner(qp.nX(), qp.nX());
			const Eigen::MatrixXd R0 = H0.bottomRightCorner(qp.nU(), qp.nU());
			const Eigen::MatrixXd S0 = H0.topRightCorner(qp.nX(), qp.nU());
			const Eigen::MatrixXd S0T = H0.bottomLeftCorner(qp.nU(), qp.nX());

			Eigen::MatrixXd A0(qp.nX(), qp.nX());
			A0 << 1, 1, 0, 1;

			Eigen::MatrixXd B0(qp.nX(), qp.nU());
			B0 << 0.5, 1.0;

			Eigen::VectorXd a0(qp.nX());
			a0 << 1, 2;

			// Stage 1
			typename Problem::StageHessianMatrix H1;
			H1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
			H1 = H1.transpose() * H1;	// Make positive definite.

			typename Problem::StageGradientVector g1;
			g1 << 0., 0., 0.;

			const Eigen::MatrixXd Q1 = H1.topLeftCorner(qp.nX(), qp.nX());
			const Eigen::MatrixXd R1 = H1.bottomRightCorner(qp.nU(), qp.nU());
			const Eigen::MatrixXd S1 = H1.topRightCorner(qp.nX(), qp.nU());
			const Eigen::MatrixXd S1T = H1.bottomLeftCorner(qp.nU(), qp.nX());

			Eigen::MatrixXd A1(qp.nX(), qp.nX());
			A1 << 1, 1, 0, 1;

			Eigen::MatrixXd B1(qp.nX(), qp.nU());
			B1 << 0.5, 1.0;

			Eigen::VectorXd a1(qp.nX());
			a1 << 1, 2;

			// Stage 2
			Eigen::MatrixXd H2(qp.nX(), qp.nX());
			H2 << 1, 2, 3, 4;
			H2 = H2.transpose() * H2;	// Make positive definite.

			typename Problem::EndStageGradientVector g2;
			g2 << 0., 0.;

			const Eigen::MatrixXd Q2 = H2.topLeftCorner(qp.nX(), qp.nX());

			// Setup QP
			set_H(qp, 0, H0);	set_g(qp, 0, g0);
			set_H(qp, 1, H1);	set_g(qp, 1, g1);
			set_Hend(qp, H2);	set_gend(qp, g2);

			typename Problem::InterStageMatrix C0;
			C0 << A0, B0;
			set_C(qp, 0, C0);
			set_c(qp, 0, a0);

			typename Problem::InterStageMatrix C1;
			C1 << A1, B1;
			set_C(qp, 1, C1);
			set_c(qp, 1, a1);
		}

		template <typename Problem>
		void problem_1(Problem& qp)
		{
			problem_0(qp);

			typename Problem::StateVector x0;
			x0 << 1., 0.;
			set_xMin(qp, 0, x0);	set_xMax(qp, 0, x0);

			set_c(qp, 0, Problem::StateVector::Zero());
			set_c(qp, 1, Problem::StateVector::Zero());
		}
	}
}
