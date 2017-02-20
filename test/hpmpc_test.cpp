/*
 * hpmpc_test.cpp
 *
 *  Created on: Jun 17, 2016
 *      Author: kotlyar
 */

#include <tmpc/qp/HPMPCProblem.hpp>
#include <tmpc/Matrix.hpp>

#include <gtest/gtest.h>

#include <vector>

TEST(HPMPCProblemTest, hpmpc_interface_works)
{
	using namespace tmpc;

	unsigned const NT = 3;
	typedef HPMPCProblem<double> Problem;
	typedef CustomMatrix<double const, unaligned, unpadded, rowMajor> MatrixMap;
	typedef CustomVector<double const, unaligned, unpadded, columnVector> VectorMap;

	std::vector<QpSize> sz;
	sz.reserve(NT);

	for (size_t k = 0; k < NT; ++k)
	{
		Rand<size_t> rand;
		sz.push_back(QpSize(rand.generate(0, 10), rand.generate(0, 10), rand.generate(0, 10)));
	}

	Problem qp(sz.begin(), sz.end());

	Rand<DynamicVector<double>> rand_vector;
	Rand<DynamicMatrix<double>> rand_matrix;

	for (size_t k = 0; k < NT; ++k)
	{
		auto const NX = sz[k].nx();
		auto const NU = sz[k].nu();
		auto const NC = sz[k].nc();
		auto const NX1 = k + 1 < NT ? sz[k + 1].nx() : 0;

		DynamicVector<double> lbx = rand_vector.generate(NX);
		DynamicVector<double> lbu = rand_vector.generate(NU);
		DynamicVector<double> ubx = rand_vector.generate(NX);
		DynamicVector<double> ubu = rand_vector.generate(NU);

		DynamicMatrix<double> H = rand_matrix.generate(NX + NC, NX + NC);
		H = trans(H) * H;	// Make positive definite.

		auto Q = submatrix(H, 0, 0, NX, NX);
		auto R = submatrix(H, NX, NX, NU, NU);
		auto S = submatrix(H, 0, NX, NX, NU);

		DynamicVector<double> q = rand_vector.generate(NX);
		DynamicVector<double> r = rand_vector.generate(NU);

		DynamicMatrix<double> A = rand_matrix.generate(NX1, NX);
		DynamicMatrix<double> B = rand_matrix.generate(NX1, NX);
		DynamicVector<double> b = rand_vector.generate(NX1);

		DynamicMatrix<double> C = rand_matrix.generate(NC, NX);
		DynamicMatrix<double> D = rand_matrix.generate(NC, NX);
		DynamicVector<double> lbd = rand_vector.generate(NC);
		DynamicVector<double> ubd = rand_vector.generate(NC);

		qp[k].set_Q(Q);
		qp[k].set_R(R);
		qp[k].set_S(S);
		qp[k].set_q(q);
		qp[k].set_r(r);

		qp[k].set_A(A);
		qp[k].set_B(B);
		qp[k].set_b(b);

		qp[k].set_C(C);
		qp[k].set_D(D);
		qp[k].set_lbd(lbd);
		qp[k].set_ubd(ubd);

		qp[k].set_lbx(lbx);
		qp[k].set_lbu(lbu);
		qp[k].set_ubx(ubx);
		qp[k].set_ubu(ubu);

		EXPECT_EQ(MatrixMap(qp.Q_data()[k], NX, NX), Q);
		EXPECT_EQ(MatrixMap(qp.R_data()[k], NU, NU), R);
		EXPECT_EQ(MatrixMap(qp.S_data()[k], NX, NU), S);
		EXPECT_EQ(VectorMap(qp.q_data()[k] + NU, NX), q);
		EXPECT_EQ(VectorMap(qp.r_data()[k] + NU, NX), r);

		EXPECT_EQ(MatrixMap(qp.A_data()[k], NX1, NX), A);
		EXPECT_EQ(MatrixMap(qp.B_data()[k], NX1, NU), B);
		EXPECT_EQ(VectorMap(qp.b_data()[k], NX1), b);

		EXPECT_EQ(MatrixMap(qp.C_data()[k], NC, NX), C);
		EXPECT_EQ(MatrixMap(qp.D_data()[k], NC, NU), D);
		EXPECT_EQ(VectorMap(qp.lg_data()[k], NC), lbd);
		EXPECT_EQ(VectorMap(qp.ug_data()[k], NC), ubd);

		EXPECT_EQ(VectorMap(qp.lb_data()[k] + NU, NX), lbx);
		EXPECT_EQ(VectorMap(qp.lb_data()[k]     , NU), lbu);
		EXPECT_EQ(VectorMap(qp.ub_data()[k] + NU, NX), ubx);
		EXPECT_EQ(VectorMap(qp.ub_data()[k]     , NU), ubu);
	}
}

