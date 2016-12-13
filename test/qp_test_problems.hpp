#pragma once

#include <stdexcept>

namespace tmpc_test
{
	namespace qp_problems
	{
		template <typename QP>
		void problem_0(QP& qp)
		{
			using namespace tmpc;

			unsigned const NX = 2;
			unsigned const NU = 1;
			unsigned const NZ = NX + NU;
			unsigned const NC = 0;
			unsigned const NCT = 0;
			unsigned const NT = 2;

			typedef Eigen::Matrix<double, NZ, NZ> StageHessianMatrix;
			typedef Eigen::Matrix<double, NZ,  1> StageGradientVector;

			if (qp.nT() != NT)
				throw std::logic_error("Invalid size of the QP problem");

			set_xu_min(qp, 0, -1.);	set_xu_max(qp, 0, 1.);
			set_xu_min(qp, 1, -1.);	set_xu_max(qp, 1, 1.);
			set_x_end_min(qp, -1.);	set_x_end_max(qp, 1.);

			// Stage 0
			StageHessianMatrix H0;
			H0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
			H0 = H0.transpose() * H0;	// Make positive definite.

			StageGradientVector g0;
			g0 << 0., 0., 0.;

			const Eigen::MatrixXd Q0 = H0.topLeftCorner(NX, NX);
			const Eigen::MatrixXd R0 = H0.bottomRightCorner(NU, NU);
			const Eigen::MatrixXd S0 = H0.topRightCorner(NX, NU);
			const Eigen::MatrixXd S0T = H0.bottomLeftCorner(NU, NX);

			Eigen::MatrixXd A0(NX, NX);
			A0 << 1, 1, 0, 1;

			Eigen::MatrixXd B0(NX, NU);
			B0 << 0.5, 1.0;

			Eigen::VectorXd a0(NX);
			a0 << 1, 2;

			// Stage 1
			StageHessianMatrix H1;
			H1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
			H1 = H1.transpose() * H1;	// Make positive definite.

			StageGradientVector g1;
			g1 << 0., 0., 0.;

			const Eigen::MatrixXd Q1 = H1.topLeftCorner(NX, NX);
			const Eigen::MatrixXd R1 = H1.bottomRightCorner(NU, NU);
			const Eigen::MatrixXd S1 = H1.topRightCorner(NX, NU);
			const Eigen::MatrixXd S1T = H1.bottomLeftCorner(NU, NX);

			Eigen::MatrixXd A1(NX, NX);
			A1 << 1, 1, 0, 1;

			Eigen::MatrixXd B1(NX, NU);
			B1 << 0.5, 1.0;

			Eigen::VectorXd a1(NX);
			a1 << 1, 2;

			// Stage 2
			Eigen::MatrixXd H2(NX, NX);
			H2 << 1, 2, 3, 4;
			H2 = H2.transpose() * H2;	// Make positive definite.

			Eigen::Matrix<double, NX, 1> g2;
			g2 << 0., 0.;

			const Eigen::MatrixXd Q2 = H2.topLeftCorner(NX, NX);

			// Setup QP
			set_H(qp, 0, H0);	set_g(qp, 0, g0);
			set_H(qp, 1, H1);	set_g(qp, 1, g1);
			qp.set_Q( 2, H2);	qp.set_q( 2, g2);

			qp.set_A(0, A0);	qp.set_B(0, B0);		qp.set_b(0, a0);
			qp.set_A(1, A1);	qp.set_B(1, B1);		qp.set_b(1, a1);
		}

		template <typename Problem>
		void problem_1(Problem& qp)
		{
			problem_0(qp);

			Eigen::Vector2d x0;
			x0 << 1., 0.;
			qp.set_x_min(0, x0);	qp.set_x_max(0, x0);

			qp.set_b(0, Eigen::Vector2d::Zero());
			qp.set_b(1, Eigen::Vector2d::Zero());
		}
	}
}
