#include <tmpc/qp/Condensing.hpp>
#include <tmpc/qp/QpOasesProblem.hpp>
#include <tmpc/kernel/eigen.hpp>
#include <tmpc/core/problem_specific.hpp>
#include <tmpc/core/RealtimeIteration.hpp>

#include "qp_test_problems.hpp"
#include "gtest_tools_eigen.hpp"

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
