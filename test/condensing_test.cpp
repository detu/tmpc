#include <tmpc/qp/Condensing.hpp>
#include <tmpc/qp/QpOasesProblem.hpp>
#include <tmpc/qp/QuadraticProblem.hpp>
#include <tmpc/Matrix.hpp>
#include <tmpc/core/problem_specific.hpp>
#include <tmpc/core/RealtimeIteration.hpp>

#include "qp_test_problems.hpp"
#include "gtest_tools_eigen.hpp"

#include <gtest/gtest.h>

#include <iostream>

using namespace tmpc;

TEST(QpSizeTest, test_CondensedQpSize)
{
	EXPECT_EQ(tmpc::condensedQpSize({
		tmpc::QpSize(2, 1, 3),
		tmpc::QpSize(4, 5, 6),
		tmpc::QpSize(2, 1, 1)
		}),
		tmpc::QpSize(2, 1 + 5 + 1, (3 + 6 + 1) + (4 + 2)));
}

TEST(QpSizeTest, test_CondensedQpSize_iteratorRange)
{
	std::array<tmpc::QpSize, 3> sz = {
		tmpc::QpSize(2, 1, 3),
		tmpc::QpSize(4, 5, 6),
		tmpc::QpSize(2, 1, 1)
		};

	EXPECT_EQ(tmpc::condensedQpSize(sz.begin(), sz.end()),
		tmpc::QpSize(2, 1 + 5 + 1, (3 + 6 + 1) + (4 + 2)));
}

class CondensingTestBlaze : public ::testing::Test
{
public:

protected:
	// Define dimensions
	struct Dimensions
	{
		static unsigned constexpr NX = 2;
		static unsigned constexpr NU = 1;
		static unsigned constexpr NW = 0;
		static unsigned constexpr NY = 0;
		static unsigned constexpr NP = 0;
		static unsigned constexpr NC = 0;
		static unsigned constexpr NCT = 0;
	};

	static unsigned constexpr NX = Dimensions::NX;
	static unsigned constexpr NU = Dimensions::NU;
	static unsigned constexpr NW = Dimensions::NW;
	static unsigned constexpr NY = Dimensions::NY;
	static unsigned constexpr NP = Dimensions::NP;
	static unsigned constexpr NC = Dimensions::NC;
	static unsigned constexpr NCT = Dimensions::NCT;
	static auto constexpr NT = 2u;
	static auto constexpr NZ = NX + NU;

	typedef double Scalar;
	typedef QuadraticProblem<Scalar> Problem;
	typedef blaze::StaticVector<Scalar, NX> StateVector;
	typedef blaze::StaticVector<Scalar, NU> InputVector;
	typedef blaze::StaticMatrix<Scalar, NX, NX> StateStateMatrix;
	typedef blaze::StaticMatrix<Scalar, NX, NU> StateInputMatrix;

	CondensingTestBlaze()
	:	size_(tmpc::RtiQpSize(NT, NX, NU, NC, NCT))
	,	qp(size_.begin(), size_.end())
	{
	}

	void SetUp() override
	{
		qp[0].set_lbx(StateVector(-1.));	qp[0].set_ubx(StateVector(1.));
		qp[0].set_lbu(InputVector(-1.));	qp[0].set_ubu(InputVector(1.));

		qp[1].set_lbx(StateVector(-1.));	qp[1].set_ubx(StateVector(1.));
		qp[1].set_lbu(InputVector(-1.));	qp[1].set_ubu(InputVector(1.));

		qp[2].set_lbx(StateVector(-1.));	qp[2].set_ubx(StateVector(1.));

		// Stage 0
		blaze::StaticMatrix<Scalar, NZ, NZ> H0 {
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9}
		};

		H0 = trans(H0) * H0;	// Make positive definite.

		const auto Q0  = submatrix(H0,  0,  0, NX, NX);
		const auto R0  = submatrix(H0, NX, NX, NU, NU);
		const auto S0  = submatrix(H0,  0, NX, NX, NU);
		const auto S0T = submatrix(H0, NX,  0, NU, NX);

		StateStateMatrix A0 {
			{1, 1},
			{0, 1}
		};

		StateInputMatrix B0 {
			{0.5},
			{1.0}
		};

		StateVector a0 {1, 2};

		// Stage 1
		blaze::StaticMatrix<Scalar, NZ, NZ> H1 {
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9}
		};
		H1 = trans(H1) * H1;	// Make positive definite.

		const auto Q1  = submatrix(H1,  0,  0, NX, NX);
		const auto R1  = submatrix(H1, NX, NX, NU, NU);
		const auto S1  = submatrix(H1,  0, NX, NX, NU);
		const auto S1T = submatrix(H1, NX,  0, NU, NX);

		StateStateMatrix A1 {
			{1, 1},
			{0, 1}
		};

		StateInputMatrix B1 {
			{0.5},
			{1.0}
		};

		StateVector a1 {1, 2};

		// Stage 2
		StateStateMatrix H2 {
			{1, 2},
			{3, 4}
		};
		H2 = trans(H2) * H2;	// Make positive definite.

		const auto Q2 = submatrix(H2, 0, 0, NX, NX);

		// Setup QP
		qp[0].set_Q(Q0);	qp[0].set_R(R0);	qp[0].set_S(S0);	qp[0].set_q(StateVector(0.));	qp[0].set_r(InputVector(0.));
		qp[1].set_Q(Q1);	qp[1].set_R(R1);	qp[1].set_S(S1);	qp[1].set_q(StateVector(0.));	qp[1].set_r(InputVector(0.));
		qp[2].set_Q(Q2);											qp[2].set_q(StateVector(0.));

		qp[0].set_A(A0);	qp[0].set_B(B0);	qp[0].set_b(a0);
		qp[1].set_A(A1);	qp[1].set_B(B1);	qp[1].set_b(a1);

		Qc_expected = trans(A0) * Q1 * A0 + trans(A0) * trans(A1) * Q2 * A1 * A0 + Q0;
		submatrix(Rc_expected, 0 * NU, 0 * NU, NU, NU) = trans(B0) * Q1 * B0 + trans(B0) * trans(A1) * Q2 * A1 * B0 + R0;
		submatrix(Rc_expected, 0 * NU, 1 * NU, NU, NU) = trans(B0) * S1 + trans(B0) * trans(A1) * Q2 * B1;
		submatrix(Rc_expected, 1 * NU, 0 * NU, NU, NU) = trans(S1) * B0 + trans(B1) * Q2 * A1 * B0;
		submatrix(Rc_expected, 1 * NU, 1 * NU, NU, NU) = trans(B1) * Q2 * B1 + R1;
		submatrix(Sc_expected, 0, 0 * NU, NX, NU) = trans(A0) * Q1 * B0 + trans(A0) * trans(A1) * Q2 * A1 * B0 + S0;
		submatrix(Sc_expected, 0, 1 * NU, NX, NU) = trans(A0) * S1 + trans(A0) * trans(A1) * Q2 * B1;

		qc_expected                        = trans(A0) * Q1 * a0				+ trans(A0) * trans(A1) * Q2 * a1		+ trans(A0) * trans(A1) * Q2 * A1 * a0;
		subvector(rc_expected, 0 * NU, NU) = trans(B0) * trans(A1) * Q2 * a1	+ trans(B0) * trans(A1) * Q2 * A1 * a0	+ trans(B0) * Q1 * a0;
		subvector(rc_expected, 1 * NU, NU) = trans(B1) * Q2 * A1 * a0			+ trans(B1) * Q2 * a1					+ trans(S1) * a0;
	}

	std::vector<tmpc::QpSize> size_;
	Problem qp;
	blaze::StaticMatrix<Scalar, NX, NX> Qc_expected;
	blaze::StaticMatrix<Scalar, NT * NU, NT * NU> Rc_expected;
	blaze::StaticMatrix<Scalar, NX, NT * NU> Sc_expected;
	blaze::StaticVector<Scalar, NX> qc_expected;
	blaze::StaticVector<Scalar, NT * NU> rc_expected;
};

unsigned constexpr CondensingTestBlaze::NX;
unsigned constexpr CondensingTestBlaze::NU;

TEST_F(CondensingTestBlaze, new_condensing_test)
{
	// Condense
	QuadraticProblem<Scalar> condensed({condensedQpSize(qpSizeIterator(qp.begin()), qpSizeIterator(qp.end()))});
	condensed.front() = condense<Scalar>(qp.begin(), qp.end());

	EXPECT_EQ(print_wrap(condensed[0].get_Q()), print_wrap(Qc_expected));
	EXPECT_EQ(print_wrap(condensed[0].get_R()), print_wrap(Rc_expected));
	EXPECT_EQ(print_wrap(condensed[0].get_S()), print_wrap(Sc_expected));
	EXPECT_EQ(print_wrap(condensed[0].get_q()), print_wrap(qc_expected));
	EXPECT_EQ(print_wrap(condensed[0].get_r()), print_wrap(rc_expected));
}
