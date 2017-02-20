#pragma once

#include <tmpc/qp/QpSize.hpp>
#include <tmpc/Matrix.hpp>

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

			typedef StaticMatrix<double, NZ, NZ> StageHessianMatrix;
			typedef StaticVector<double, NZ> StageGradientVector;

			if (qp.nT() != NT)
				throw std::logic_error("Invalid size of the QP problem");

			set_xu_min(qp, 0, -1.);	set_xu_max(qp, 0, 1.);
			set_xu_min(qp, 1, -1.);	set_xu_max(qp, 1, 1.);
			set_x_end_min(qp, -1.);	set_x_end_max(qp, 1.);

			// Stage 0
			StageHessianMatrix H0 {
					{1, 2, 3},
					{4, 5, 6},
					{7, 8, 9}
			};
			H0 = trans(H0) * H0;	// Make positive definite.

			StageGradientVector g0 {0., 0., 0.};

			const DynamicMatrix<double> Q0 = submatrix(H0, 0, 0, NX, NX);
			const DynamicMatrix<double> R0 = submatrix(H0, NX, NX, NU, NU);
			const DynamicMatrix<double> S0 = submatrix(H0, 0, NX, NX, NU);
			const DynamicMatrix<double> S0T = submatrix(H0, NX, 0, NU, NX);

			DynamicMatrix<double> const A0 {{1, 1}, {0, 1}};
			DynamicMatrix<double> const B0 {{0.5}, {1.0}};
			DynamicVector<double> a0 {1, 2};

			// Stage 1
			StageHessianMatrix H1 {
					{1, 2, 3},
					{4, 5, 6},
					{7, 8, 9}
			};
			H1 = trans(H1) * H1;	// Make positive definite.

			StageGradientVector g1 {0., 0., 0.};

			const DynamicMatrix<double> Q1 = submatrix(H1, 0, 0, NX, NX);
			const DynamicMatrix<double> R1 = submatrix(H1, NX, NX, NU, NU);
			const DynamicMatrix<double> S1 = submatrix(H1, 0, NX, NX, NU);
			const DynamicMatrix<double> S1T = submatrix(H1, NX, 0, NU, NX);

			DynamicMatrix<double> const A1 {{1, 1}, {0, 1}};
			DynamicMatrix<double> const B1 {{0.5}, {1.0}};
			DynamicVector<double> const a1 {1, 2};

			// Stage 2
			DynamicMatrix<double> H2 {{1, 2}, {3, 4}};
			H2 = trans(H2) * H2;	// Make positive definite.

			StaticVector<double, NX> const g2 {0., 0.};

			const DynamicMatrix<double> Q2 = submatrix(H2, 0, 0, NX, NX);

			// Setup QP
			set_H(qp, 0, H0);	set_g(qp, 0, g0);
			set_H(qp, 1, H1);	set_g(qp, 1, g1);
			qp.set_Q( 2, H2);	qp.set_q( 2, g2);

			qp.set_A(0, A0);	qp.set_B(0, B0);		qp.set_b(0, a0);
			qp.set_A(1, A1);	qp.set_B(1, B1);		qp.set_b(1, a1);
		}

		template <typename QP>
		QP problem_0()
		{
			using namespace tmpc;

			unsigned const NX = 2;
			unsigned const NU = 1;
			unsigned const NC = 0;
			unsigned const NCT = 0;
			unsigned const NT = 2;

			QP qp({QpSize(NX, NU, NC), QpSize(NX, NU, NC), QpSize(NX, 0, NCT)});
			problem_0(qp);

			return qp;
		}

		template <typename Problem>
		void problem_1(Problem& qp)
		{
			problem_0(qp);

			tmpc::StaticVector<double, 2> x0 {1., 0.};
			qp.set_x_min(0, x0);	qp.set_x_max(0, x0);

			qp.set_b(0, tmpc::StaticVector<double, 2>(0.));
			qp.set_b(1, tmpc::StaticVector<double, 2>(0.));
		}
	}
}
