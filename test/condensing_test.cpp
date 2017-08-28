#include <tmpc/qp/Condensing.hpp>
#include <tmpc/qp/QpOasesProblem.hpp>
#include <tmpc/qp/QuadraticProblem.hpp>
#include <tmpc/kernel/eigen.hpp>
#include <tmpc/core/problem_specific.hpp>
#include <tmpc/core/RealtimeIteration.hpp>
#include <tmpc/test_tools.hpp>

#include "qp_test_problems.hpp"

#include <gtest/gtest.h>

#include <iostream>

class CondensingTest : public ::testing::Test
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

	typedef tmpc::EigenKernel<double> KBase;

	// Define a kernel augmented with problem-specific types
	struct K : KBase, tmpc::ProblemSpecific<KBase, Dimensions> {};

	static unsigned constexpr NX = Dimensions::NX;
	static unsigned constexpr NU = Dimensions::NU;
	static unsigned constexpr NW = Dimensions::NW;
	static unsigned constexpr NY = Dimensions::NY;
	static unsigned constexpr NP = Dimensions::NP;
	static unsigned constexpr NC = Dimensions::NC;
	static unsigned constexpr NCT = Dimensions::NCT;
	static auto constexpr NT = 2u;
	static auto constexpr NZ = K::NX + K::NU;

	typedef tmpc::QpOasesProblem Problem;

	CondensingTest()
	:	qp(tmpc::RtiQpSize(NT, NX, NU, NC, NCT))
	,	Hc_expected(NX + NT * NU, NX + NT * NU)
	,	gc_expected(NX + NT * NU)
	{
	}

	void SetUp() override
	{
		qp.set_x_min(0, K::StateVector::Constant(-1.));	qp.set_x_max(0, K::StateVector::Constant(1.));
		qp.set_u_min(0, K::InputVector::Constant(-1.));	qp.set_u_max(0, K::InputVector::Constant(1.));

		qp.set_x_min(1, K::StateVector::Constant(-1.));	qp.set_x_max(1, K::StateVector::Constant(1.));
		qp.set_u_min(1, K::InputVector::Constant(-1.));	qp.set_u_max(1, K::InputVector::Constant(1.));

		qp.set_x_min(2, K::StateVector::Constant(-1.));	qp.set_x_max(2, K::StateVector::Constant(1.));

		// Stage 0
		K::Matrix<NZ, NZ> H0;
		H0 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		H0 = H0.transpose() * H0;	// Make positive definite.

		const Eigen::MatrixXd Q0 = H0.topLeftCorner(NX, NX);
		const Eigen::MatrixXd R0 = H0.bottomRightCorner(NU, NU);
		const Eigen::MatrixXd S0 = H0.topRightCorner(NX, NU);
		const Eigen::MatrixXd S0T = H0.bottomLeftCorner(NU, NX);

		K::StateStateMatrix A0;
		A0 << 1, 1, 0, 1;

		K::StateInputMatrix B0(NX, NU);
		B0 << 0.5, 1.0;

		Eigen::VectorXd a0(NX);
		a0 << 1, 2;

		// Stage 1
		K::Matrix<NZ, NZ> H1;
		H1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
		H1 = H1.transpose() * H1;	// Make positive definite.

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
		K::StateStateMatrix H2;
		H2 << 1, 2, 3, 4;
		H2 = H2.transpose() * H2;	// Make positive definite.

		const Eigen::MatrixXd Q2 = H2.topLeftCorner(NX, NX);

		// Setup QP
		set_H(qp, 0, H0);	qp.set_q(0, K::StateVector::Zero());	qp.set_r(0, K::InputVector::Zero());
		set_H(qp, 1, H1);	qp.set_q(1, K::StateVector::Zero());	qp.set_r(1, K::InputVector::Zero());
		set_Q_end(qp, H2);	qp.set_q(2, K::StateVector::Zero());

		qp.set_A(0, A0);	qp.set_B(0, B0);	qp.set_b(0, a0);
		qp.set_A(1, A1);	qp.set_B(1, B1);	qp.set_b(1, a1);

		Hc_expected <<
				A0.transpose() * Q1 * A0 + A0.transpose() * A1.transpose() * Q2 * A1 * A0 + Q0,				A0.transpose() * Q1 * B0 + A0.transpose() * A1.transpose() * Q2 * A1 * B0 + S0,	A0.transpose() * S1 + A0.transpose() * A1.transpose() * Q2 * B1,
				B0.transpose() * Q1 * A0 + B0.transpose() * A1.transpose() * Q2 * A1 * A0 + S0.transpose(), B0.transpose() * Q1 * B0 + B0.transpose() * A1.transpose() * Q2 * A1 * B0 + R0, B0.transpose() * S1 + B0.transpose() * A1.transpose() * Q2 * B1,
				S1.transpose() * A0 + B1.transpose() * Q2 * A1 * A0,										S1.transpose() * B0 + B1.transpose() * Q2 * A1 * B0,							B1.transpose() * Q2 * B1 + R1;

		gc_expected <<
				A0.transpose() * Q1 * a0					+ A0.transpose() * A1.transpose() * Q2 * a1			+ A0.transpose() * A1.transpose() * Q2 * A1 * a0,
				B0.transpose() * A1.transpose() * Q2 * a1	+ B0.transpose() * A1.transpose() * Q2 * A1 * a0	+ B0.transpose() * Q1 * a0,
				B1.transpose() * Q2 * A1 * a0				+ B1.transpose() * Q2 * a1							+ S1.transpose() * a0;
	}

	Problem qp;
	Eigen::MatrixXd Hc_expected;
	Eigen::VectorXd gc_expected;
};

unsigned constexpr CondensingTest::NX;
unsigned constexpr CondensingTest::NU;

/*
TEST_F(CondensingTest, condensing_test)
{
	// Condense
	tmpc::QpOasesProblem condensed({tmpc::condensedQpSize(tmpc::qpSizeIterator(qp.begin()), tmpc::qpSizeIterator(qp.end()))});
	tmpc::Condense<K, Dimensions>(qp, condensed);

	EXPECT_EQ(print_wrap(condensed.H()), print_wrap(Hc_expected));
	EXPECT_EQ(print_wrap(condensed.g()), print_wrap(gc_expected));
}
*/

TEST_F(CondensingTest, new_condensing_test)
{
	// Condense
	tmpc::QpOasesProblem condensed({tmpc::condensedQpSize(tmpc::qpSizeIterator(qp.begin()), tmpc::qpSizeIterator(qp.end()))});
	tmpc::Condense<K>(qp.begin(), qp.end(), condensed.front());

	EXPECT_EQ(print_wrap(condensed.H()), print_wrap(Hc_expected));
	EXPECT_EQ(print_wrap(condensed.g()), print_wrap(gc_expected));
}

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
	typedef tmpc::QuadraticProblemBlaze<Scalar> Problem;
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
		qp[0].lbx() = StateVector(-1.);	qp[0].ubx() = StateVector(1.);
		qp[0].lbu() = InputVector(-1.);	qp[0].ubu() = InputVector(1.);

		qp[1].lbx() = StateVector(-1.);	qp[1].ubx() = StateVector(1.);
		qp[1].lbu() = InputVector(-1.);	qp[1].ubu() = InputVector(1.);

		qp[2].lbx() = StateVector(-1.);	qp[2].ubx() = StateVector(1.);

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
		qp[0].Q() = (Q0);	qp[0].R() = (R0);	qp[0].S() = (S0);	qp[0].q() = (StateVector(0.));	qp[0].r() = (InputVector(0.));
		qp[1].Q() = (Q1);	qp[1].R() = (R1);	qp[1].S() = (S1);	qp[1].q() = (StateVector(0.));	qp[1].r() = (InputVector(0.));
		qp[2].Q() = (Q2);											qp[2].q() = (StateVector(0.));

		qp[0].A() = (A0);	qp[0].B() = (B0);	qp[0].b() = (a0);
		qp[1].A() = (A1);	qp[1].B() = (B1);	qp[1].b() = (a1);

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
	tmpc::QuadraticProblemBlaze<Scalar> condensed({tmpc::condensedQpSize(tmpc::qpSizeIterator(qp.begin()), tmpc::qpSizeIterator(qp.end()))});
	condensed.front() = tmpc::condenseBlaze<Scalar>(qp.begin(), qp.end());

	EXPECT_EQ(print_wrap(condensed[0].Q()), print_wrap(Qc_expected));
	EXPECT_EQ(print_wrap(condensed[0].R()), print_wrap(Rc_expected));
	EXPECT_EQ(print_wrap(condensed[0].S()), print_wrap(Sc_expected));
	EXPECT_EQ(print_wrap(condensed[0].q()), print_wrap(qc_expected));
	EXPECT_EQ(print_wrap(condensed[0].r()), print_wrap(rc_expected));
}
