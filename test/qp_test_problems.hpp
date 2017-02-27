#pragma once

#include <tmpc/qp/QpSize.hpp>
#include <tmpc/Matrix.hpp>

#include <stdexcept>

namespace tmpc_test
{
	namespace qp_problems
	{
		template <typename QP>
		QP problem_0()
		{
			using namespace tmpc;
			
			unsigned const NX = 2;
			unsigned const NU = 1;
			unsigned const NZ = NX + NU;
			unsigned const NC = 0;
			unsigned const NCT = 0;
			unsigned const NT = 2;

			typedef StaticMatrix<double, NZ, NZ> StageHessianMatrix;

			const auto sz = RtiQpSize(NT, NX, NU, NC, NCT);
			QP qp(sz.begin(), sz.end());
			
			qp[0].set_lbx(-1.);	qp[0].set_lbu(-1.);	qp[0].set_ubx(1.);	qp[0].set_ubu(1.);
			qp[1].set_lbx(-1.);	qp[1].set_lbu(-1.);	qp[1].set_ubx(1.);	qp[1].set_ubu(1.);
			qp[2].set_lbx(-1.);						qp[2].set_ubx(1.);

			// Stage 0
			StageHessianMatrix H0 {
					{1, 2, 3},
					{4, 5, 6},
					{7, 8, 9}
			};
			H0 = trans(H0) * H0;	// Make positive definite.

			StaticVector<double, NX> const q0 {0., 0.};
			StaticVector<double, NU> const r0 {0.};

			const DynamicMatrix<double> Q0 = submatrix(H0, 0, 0, NX, NX);
			const DynamicMatrix<double> R0 = submatrix(H0, NX, NX, NU, NU);
			const DynamicMatrix<double> S0 = submatrix(H0, 0, NX, NX, NU);
			const DynamicMatrix<double> S0T = submatrix(H0, NX, 0, NU, NX);

			DynamicMatrix<double> const A0 {{1., 1.}, {0., 1.}};
			DynamicMatrix<double> const B0 {{0.5}, {1.0}};
			DynamicVector<double> a0 {1., 2.};

			// Stage 1
			StageHessianMatrix H1 {
					{1., 2., 3.},
					{4., 5., 6.},
					{7., 8., 9.}
			};
			H1 = trans(H1) * H1;	// Make positive definite.

			StaticVector<double, NX> const q1 {0., 0.};
			StaticVector<double, NU> const r1 {0.};

			const DynamicMatrix<double> Q1 = submatrix(H1, 0, 0, NX, NX);
			const DynamicMatrix<double> R1 = submatrix(H1, NX, NX, NU, NU);
			const DynamicMatrix<double> S1 = submatrix(H1, 0, NX, NX, NU);
			const DynamicMatrix<double> S1T = submatrix(H1, NX, 0, NU, NX);

			DynamicMatrix<double> const A1 {{1., 1.}, {0., 1.}};
			DynamicMatrix<double> const B1 {{0.5}, {1.0}};
			DynamicVector<double> const a1 {1., 2.};

			// Stage 2
			DynamicMatrix<double> H2 {{1., 2.}, {3., 4.}};
			H2 = trans(H2) * H2;	// Make positive definite.

			StaticVector<double, NX> const q2 {0., 0.};

			const DynamicMatrix<double> Q2 = submatrix(H2, 0, 0, NX, NX);

			// Setup QP
			qp[0].set_Q(Q0);	qp[0].set_R(R0);	qp[0].set_S(S0);	qp[0].set_q(q0);	qp[0].set_r(r0);
			qp[1].set_Q(Q1);	qp[1].set_R(R1);	qp[1].set_S(S1);	qp[1].set_q(q1);	qp[1].set_r(r1);
			qp[2].set_Q(Q2);											qp[2].set_q(q2);

			qp[0].set_A(A0);	
			qp[0].set_B(B0);		
			qp[0].set_b(a0);
			
			qp[1].set_A(A1);	
			qp[1].set_B(B1);		
			qp[1].set_b(a1);

			return std::move(qp);
		}

		template <typename Problem>
		Problem problem_1()
		{
			Problem qp = problem_0<Problem>();

			tmpc::StaticVector<double, 2> x0 {1., 0.};
			qp[0].set_lbx(x0);	qp[0].set_ubx(x0);

			qp[0].set_b(0.);
			qp[1].set_b(0.);

			return std::move(qp);
		}
	}
}
