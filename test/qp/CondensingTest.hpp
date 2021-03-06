#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	template <typename CA>
	class CondensingTest 
	: 	public Test
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

		using CondensingAlgorithm = CA;
		using Kernel = typename CondensingAlgorithm::Kernel;
		using Real = typename Kernel::Real;
		using StateVector = StaticVector<Kernel, NX>;
		using InputVector = StaticVector<Kernel, NU>;
		using StateStateMatrix = StaticMatrix<Kernel, NX, NX>;
		using StateInputMatrix = StaticMatrix<Kernel, NX, NU>;

		CondensingTest()
		:	size_ {mpcOcpSize(NT, NX, NU, NC, NCT)}
		{
			for (size_t k = 0; k < size_.size(); ++k)
				qp.emplace_back(size_[k], k + 1 < size_.size() ? size_[k + 1].nx() : 0);
		}

		void SetUp() override
		{
			qp[0].lbx(StateVector(-1.));	qp[0].ubx(StateVector(1.));
			qp[0].lbu(InputVector(-1.));	qp[0].ubu(InputVector(1.));

			qp[1].lbx(StateVector(-1.));	qp[1].ubx(StateVector(1.));
			qp[1].lbu(InputVector(-1.));	qp[1].ubu(InputVector(1.));

			qp[2].lbx(StateVector(-1.));	qp[2].ubx(StateVector(1.));

			// Stage 0
			StaticMatrix<Kernel, NZ, NZ> H0 {
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9}
			};

			H0 = trans(H0) * H0;	// Make positive definite.

			const auto Q0  = submatrix(H0,  0,  0, NX, NX);
			const auto R0  = submatrix(H0, NX, NX, NU, NU);
			const auto S0  = submatrix(H0,  0, NX, NX, NU);

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
			StaticMatrix<Kernel, NZ, NZ> H1 {
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9}
			};
			H1 = trans(H1) * H1;	// Make positive definite.

			const auto Q1  = submatrix(H1,  0,  0, NX, NX);
			const auto R1  = submatrix(H1, NX, NX, NU, NU);
			const auto S1  = submatrix(H1,  0, NX, NX, NU);

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

		std::vector<DynamicOcpSize> size_;
		std::vector<OcpQpStage<Kernel>> qp;
		StaticMatrix<Kernel, NX, NX> Qc_expected;
		StaticMatrix<Kernel, NT * NU, NT * NU> Rc_expected;
		StaticMatrix<Kernel, NX, NT * NU> Sc_expected;
		StaticVector<Kernel, NX> qc_expected;
		StaticVector<Kernel, NT * NU> rc_expected;
		StaticVector<Kernel, NX> lbx_expected;
		StaticVector<Kernel, NX> ubx_expected;
		StaticVector<Kernel, NT * NU> lbu_expected;
		StaticVector<Kernel, NT * NU> ubu_expected;
	};


	TYPED_TEST_SUITE_P(CondensingTest);

	
	TYPED_TEST_P(CondensingTest, testCondensing)
	{
		using Kernel = typename TestFixture::Kernel;

		// Condense
		typename TestFixture::CondensingAlgorithm condensing(sizeBegin(this->qp), sizeEnd(this->qp));
		OcpQpStage<Kernel> const condensed = condensing(this->qp.begin(), this->qp.end());

		EXPECT_EQ(forcePrint(condensed.Q()), forcePrint(this->Qc_expected));
		EXPECT_EQ(forcePrint(condensed.R()), forcePrint(this->Rc_expected));
		EXPECT_EQ(forcePrint(condensed.S()), forcePrint(this->Sc_expected));
		EXPECT_EQ(forcePrint(condensed.q()), forcePrint(this->qc_expected));
		EXPECT_EQ(forcePrint(condensed.r()), forcePrint(this->rc_expected));
		EXPECT_EQ(forcePrint(condensed.lbx()), forcePrint(this->lbx_expected));
		EXPECT_EQ(forcePrint(condensed.ubx()), forcePrint(this->ubx_expected));
		EXPECT_EQ(forcePrint(condensed.lbu()), forcePrint(this->lbu_expected));
		EXPECT_EQ(forcePrint(condensed.ubu()), forcePrint(this->ubu_expected));
	}
	

	REGISTER_TYPED_TEST_SUITE_P(CondensingTest,
		testCondensing
	);
}