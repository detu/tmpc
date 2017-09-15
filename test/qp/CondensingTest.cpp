#include <tmpc/qp/Condensing.hpp>
#include <tmpc/qp/QuadraticProblem.hpp>
#include <tmpc/mpc/MpcQpSize.hpp>

#include <tmpc/BlazeKernel.hpp>

#include "../gtest_tools_eigen.hpp"

#include <gtest/gtest.h>

namespace tmpc :: testing
{
	TEST(QpSizeTest, testCondensedQpSize)
	{
		EXPECT_EQ(tmpc::condensedQpSize({
			tmpc::QpSize(2, 1, 3),
			tmpc::QpSize(4, 5, 6),
			tmpc::QpSize(2, 1, 1)
			}),
			tmpc::QpSize(2, 1 + 5 + 1, (3 + 6 + 1) + (4 + 2)));
	}

	TEST(QpSizeTest, testCondensedQpSizeIteratorRange)
	{
		std::array<tmpc::QpSize, 3> sz = {
			tmpc::QpSize(2, 1, 3),
			tmpc::QpSize(4, 5, 6),
			tmpc::QpSize(2, 1, 1)
			};

		EXPECT_EQ(tmpc::condensedQpSize(sz.begin(), sz.end()),
			tmpc::QpSize(2, 1 + 5 + 1, (3 + 6 + 1) + (4 + 2)));
	}

	class CondensingTest : public ::testing::Test
	{
	public:

	protected:
		// Define dimensions
		static unsigned constexpr NX = 2;
		static unsigned constexpr NU = 1;
		static unsigned constexpr NW = 0;
		static unsigned constexpr NY = 0;
		static unsigned constexpr NP = 0;
		static unsigned constexpr NC = 0;
		static unsigned constexpr NCT = 0;
		static auto constexpr NT = 2u;
		static auto constexpr NZ = NX + NU;

		using Real = double;
		using Kernel = BlazeKernel<Real>;
		using Problem = QuadraticProblem<Kernel>;
		using StateVector = Kernel::StaticVector<NX>;
		using InputVector = Kernel::StaticVector<NU>;
		using StateStateMatrix = Kernel::StaticMatrix<NX, NX>;
		using StateInputMatrix = Kernel::StaticMatrix<NX, NU>;

		CondensingTest()
		:	size_(tmpc::mpcQpSize(NT, NX, NU, NC, NCT))
		,	qp(size_.begin(), size_.end())
		{
		}

		void SetUp() override
		{
			qp[0].lbx(StateVector(-1.));	qp[0].ubx(StateVector(1.));
			qp[0].lbu(InputVector(-1.));	qp[0].ubu(InputVector(1.));

			qp[1].lbx(StateVector(-1.));	qp[1].ubx(StateVector(1.));
			qp[1].lbu(InputVector(-1.));	qp[1].ubu(InputVector(1.));

			qp[2].lbx(StateVector(-1.));	qp[2].ubx(StateVector(1.));

			// Stage 0
			Kernel::StaticMatrix<NZ, NZ> H0 {
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
			Kernel::StaticMatrix<NZ, NZ> H1 {
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
			qp[0].Q(Q0);	qp[0].R(R0);	qp[0].S(S0);	qp[0].q(StateVector(0.));	qp[0].r(InputVector(0.));
			qp[1].Q(Q1);	qp[1].R(R1);	qp[1].S(S1);	qp[1].q(StateVector(0.));	qp[1].r(InputVector(0.));
			qp[2].Q(Q2);											qp[2].q(StateVector(0.));

			qp[0].A(A0);	qp[0].B(B0);	qp[0].b(a0);
			qp[1].A(A1);	qp[1].B(B1);	qp[1].b(a1);

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

			lbx_expected = qp[0].lbx();
			ubx_expected = qp[0].ubx();
			
			subvector(lbu_expected, 0 * NU, NU) = qp[0].lbu();	subvector(ubu_expected, 0 * NU, NU) = qp[0].ubu();
			subvector(lbu_expected, 1 * NU, NU) = qp[1].lbu();	subvector(ubu_expected, 1 * NU, NU) = qp[1].ubu();
		}

		std::vector<tmpc::QpSize> size_;
		Problem qp;
		Kernel::StaticMatrix<NX, NX> Qc_expected;
		Kernel::StaticMatrix<NT * NU, NT * NU> Rc_expected;
		Kernel::StaticMatrix<NX, NT * NU> Sc_expected;
		Kernel::StaticVector<NX> qc_expected;
		Kernel::StaticVector<NT * NU> rc_expected;
		Kernel::StaticVector<NX> lbx_expected;
		Kernel::StaticVector<NX> ubx_expected;
		Kernel::StaticVector<NT * NU> lbu_expected;
		Kernel::StaticVector<NT * NU> ubu_expected;
	};

	unsigned constexpr CondensingTest::NX;
	unsigned constexpr CondensingTest::NU;

	TEST_F(CondensingTest, testCondensing)
	{
		// Condense
		Condensing<Kernel> condensing(sizeBegin(qp), sizeEnd(qp));
		QuadraticProblemStage<Kernel> const condensed = condensing(qp.begin(), qp.end());

		EXPECT_EQ(print_wrap(condensed.Q()), print_wrap(Qc_expected));
		EXPECT_EQ(print_wrap(condensed.R()), print_wrap(Rc_expected));
		EXPECT_EQ(print_wrap(condensed.S()), print_wrap(Sc_expected));
		EXPECT_EQ(print_wrap(condensed.q()), print_wrap(qc_expected));
		EXPECT_EQ(print_wrap(condensed.r()), print_wrap(rc_expected));
		EXPECT_EQ(print_wrap(condensed.lbx()), print_wrap(lbx_expected));
		EXPECT_EQ(print_wrap(condensed.ubx()), print_wrap(ubx_expected));
		EXPECT_EQ(print_wrap(condensed.lbu()), print_wrap(lbu_expected));
		EXPECT_EQ(print_wrap(condensed.ubu()), print_wrap(ubu_expected));
	}
}