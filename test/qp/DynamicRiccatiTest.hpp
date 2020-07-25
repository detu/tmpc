#pragma once

#include <tmpc/qp/DynamicOcpQp.hpp>
#include <tmpc/qp/Randomize.hpp>
#include <tmpc/qp/KktValue.hpp>
#include <tmpc/ocp/DynamicOcpSolution.hpp>
#include <tmpc/Testing.hpp>
#include <tmpc/testing/qp/IsSolution.hpp>

#include <blaze/Math.h>


namespace tmpc :: testing
{
    template <typename RiccatiImpl_>
    class DynamicRiccatiTest 
    :   public Test
    {
    public:
        using Qp = DynamicOcpQp<double>;
		using Solution = DynamicOcpSolution<double>;
		using RiccatiImpl = RiccatiImpl_;

    protected:
    };


    TYPED_TEST_SUITE_P(DynamicRiccatiTest);


    TYPED_TEST_P(DynamicRiccatiTest, test_riccati0)
	{
		unsigned const NX = 2;
		unsigned const NU = 1;
		unsigned const NZ = NX + NU;
		unsigned const NC = 0;
		unsigned const NCT = 0;
		unsigned const NT = 2;

		using blaze::StaticMatrix;
		using blaze::StaticVector;
		using blaze::DynamicVector;
		using blaze::DynamicMatrix;
		using StageHessianMatrix = StaticMatrix<double, NZ, NZ>;

		OcpTree const g(NT + 1);
		DynamicOcpSize const sz {g, NX, NU, NC, 0, NCT, false};
		typename TestFixture::Qp qp {sz};
		typename TestFixture::Solution sol {sz};

		auto const e0 = out_edges(0, g).front();
		auto const e1 = out_edges(1, g).front();
		
		qp.lx(0, DynamicVector<double>(NX, -1.));
		qp.lu(0, DynamicVector<double>(NU, -1.));
		qp.ux(0, DynamicVector<double>(NX, 1.));
		qp.uu(0, DynamicVector<double>(NU, 1.));

		qp.lx(1, DynamicVector<double>(NX, -1.));
		qp.lu(1, DynamicVector<double>(NU, -1.));
		qp.ux(1, DynamicVector<double>(NX, 1.));
		qp.uu(1, DynamicVector<double>(NU, 1.));

		qp.lx(2, DynamicVector<double>(NX, -1.));
		qp.ux(2, DynamicVector<double>(NX, 1.));
		
		// Stage 0
		StageHessianMatrix H0 {
			{67,   78,   90},
			{78,   94,  108},
			{90,  108,  127}
		};

		StaticVector<double, NX> const q0 {0., 0.};
		StaticVector<double, NU> const r0 {0.};

		const DynamicMatrix<double> Q0 = submatrix(H0, 0, 0, NX, NX);
		const DynamicMatrix<double> R0 = submatrix(H0, NX, NX, NU, NU);
		const DynamicMatrix<double> S0T = submatrix(H0, 0, NX, NX, NU);
		const DynamicMatrix<double> S0 = submatrix(H0, NX, 0, NU, NX);

		DynamicMatrix<double> const A0 {{1., 1.}, {0., 1.}};
		DynamicMatrix<double> const B0 {{0.5}, {1.0}};
		DynamicVector<double> a0 {1., 2.};

		// Stage 1
		StageHessianMatrix H1 {
			{67,   78,   90},
			{78,   94,  108},
			{90,  108,  127}
		};

		StaticVector<double, NX> const q1 {0., 0.};
		StaticVector<double, NU> const r1 {0.};

		const DynamicMatrix<double> Q1 = submatrix(H1, 0, 0, NX, NX);
		const DynamicMatrix<double> R1 = submatrix(H1, NX, NX, NU, NU);
		const DynamicMatrix<double> S1T = submatrix(H1, 0, NX, NX, NU);
		const DynamicMatrix<double> S1 = submatrix(H1, NX, 0, NU, NX);

		DynamicMatrix<double> const A1 {{1., 1.}, {0., 1.}};
		DynamicMatrix<double> const B1 {{0.5}, {1.0}};
		DynamicVector<double> const a1 {1., 2.};

		// Stage 2
		DynamicMatrix<double> H2 {{1., 2.}, {3., 4.}};
		H2 = trans(H2) * H2;	// Make positive definite.

		StaticVector<double, NX> const q2 {0., 0.};

		const DynamicMatrix<double> Q2 = submatrix(H2, 0, 0, NX, NX);

		// Setup QP
		qp.Q(0, Q0);	
		qp.R(0, R0);	
		qp.S(0, S0);	
		qp.q(0, q0);	
		qp.r(0, r0);

		qp.Q(1, Q1);	
		qp.R(1, R1);	
		qp.S(1, S1);	
		qp.q(1, q1);	
		qp.r(1, r1);

		qp.Q(2, Q2);	
		qp.q(2, q2);

		qp.A(e0, A0);	
		qp.B(e0, B0);	
		qp.b(e0, a0);

		qp.A(e1, A1);	
		qp.B(e1, B1);	
		qp.b(e1, a1);

		typename TestFixture::RiccatiImpl riccati {sz};
		riccati(qp, sol);

		EXPECT_TRUE(approxEqual(sol.x(0), (DynamicVector<double> {4.2376727217537882, -3.2166970575479454}), 1e-6));
		EXPECT_TRUE(approxEqual(sol.u(0), (DynamicVector<double> {-0.34889319983994238}), 1e-6));
		EXPECT_TRUE(approxEqual(sol.pi(e0), (DynamicVector<double> {-1.621313883169287, 11.132830578594472}), 1e-6));

		EXPECT_TRUE(approxEqual(sol.x(1), (DynamicVector<double> {1.8465290642858716, -1.5655902573878877}), 1e-6));
		EXPECT_TRUE(approxEqual(sol.u(1), (DynamicVector<double> {-0.20287931724419153}), 1e-6));
		EXPECT_TRUE(approxEqual(sol.pi(e1), (DynamicVector<double> {15.036417437909771, 21.143596583220845}), 1e-6));

		EXPECT_TRUE(approxEqual(sol.x(2), (DynamicVector<double> {1.1794991482758881, 0.23153042536792068}), 1e-6));
	}


	TYPED_TEST_P(DynamicRiccatiTest, test_riccati1)
	{
		unsigned const NX = 2;
		unsigned const NU = 1;
		unsigned const NZ = NX + NU;
		unsigned const NC = 0;
		unsigned const NCT = 0;
		unsigned const NT = 2;

		using blaze::StaticMatrix;
		using blaze::StaticVector;
		using blaze::DynamicVector;
		using blaze::DynamicMatrix;
		using StageHessianMatrix = StaticMatrix<double, NZ, NZ>;

		OcpTree const g(NT + 1);
		DynamicOcpSize const sz {g, NX, NU, NC, 0, NCT, true};
		typename TestFixture::Qp qp {sz};
		typename TestFixture::Solution sol {sz};

		auto const e0 = out_edges(0, g).front();
		auto const e1 = out_edges(1, g).front();
		
		qp.lx(0, DynamicVector<double>(0, -1.));
		qp.lu(0, DynamicVector<double>(NU, -1.));
		qp.ux(0, DynamicVector<double>(0, 1.));
		qp.uu(0, DynamicVector<double>(NU, 1.));

		qp.lx(1, DynamicVector<double>(NX, -1.));
		qp.lu(1, DynamicVector<double>(NU, -1.));
		qp.ux(1, DynamicVector<double>(NX, 1.));
		qp.uu(1, DynamicVector<double>(NU, 1.));

		qp.lx(2, DynamicVector<double>(NX, -1.));
		qp.ux(2, DynamicVector<double>(NX, 1.));
		
		// Stage 0
		DynamicMatrix<double> H0 {
			{94,  108},
			{108,  127}
		};

		DynamicVector<double> const q0 {};
		StaticVector<double, NU> const r0 {0.};

		const DynamicMatrix<double> Q0 = submatrix(H0, 0, 0, 0, 0);
		const DynamicMatrix<double> R0 = submatrix(H0, 0, 0, NU, NU);
		const DynamicMatrix<double> S0T = submatrix(H0, 0, 0, 0, NU);
		const DynamicMatrix<double> S0 = submatrix(H0, 0, 0, NU, 0);

		DynamicMatrix<double> const A0 {
			{},
			{}
		};
		DynamicMatrix<double> const B0 {{0.5}, {1.0}};
		DynamicVector<double> a0 {1., 2.};

		// Stage 1
		StageHessianMatrix H1 {
			{67,   78,   90},
			{78,   94,  108},
			{90,  108,  127}
		};

		StaticVector<double, NX> const q1 {0., 0.};
		StaticVector<double, NU> const r1 {0.};

		const DynamicMatrix<double> Q1 = submatrix(H1, 0, 0, NX, NX);
		const DynamicMatrix<double> R1 = submatrix(H1, NX, NX, NU, NU);
		const DynamicMatrix<double> S1T = submatrix(H1, 0, NX, NX, NU);
		const DynamicMatrix<double> S1 = submatrix(H1, NX, 0, NU, NX);

		DynamicMatrix<double> const A1 {{1., 1.}, {0., 1.}};
		DynamicMatrix<double> const B1 {{0.5}, {1.0}};
		DynamicVector<double> const a1 {1., 2.};

		// Stage 2
		DynamicMatrix<double> H2 {{1., 2.}, {3., 4.}};
		H2 = trans(H2) * H2;	// Make positive definite.

		StaticVector<double, NX> const q2 {0., 0.};

		const DynamicMatrix<double> Q2 = submatrix(H2, 0, 0, NX, NX);

		// Setup QP
		qp.Q(0, Q0);	
		qp.R(0, R0);	
		qp.S(0, S0);	
		qp.q(0, q0);	
		qp.r(0, r0);

		qp.Q(1, Q1);	
		qp.R(1, R1);	
		qp.S(1, S1);	
		qp.q(1, q1);	
		qp.r(1, r1);

		qp.Q(2, Q2);	
		qp.q(2, q2);

		qp.A(e0, A0);	
		qp.B(e0, B0);	
		qp.b(e0, a0);

		qp.A(e1, A1);	
		qp.B(e1, B1);	
		qp.b(e1, a1);

		typename TestFixture::RiccatiImpl riccati {sz};
		riccati(qp, sol);

		EXPECT_TRUE(approxEqual(sol.x(0), (DynamicVector<double> {}), 1e-6));
		EXPECT_TRUE(approxEqual(sol.u(0), (DynamicVector<double> {-0.32198212467473536}), 1e-6));
		//EXPECT_TRUE(approxEqual(sol.pi(e0), (DynamicVector<double> {-8.1496398536788472, 34.341139646264558}), 1e-6));

		EXPECT_TRUE(approxEqual(sol.x(1), (DynamicVector<double> {0.83900893766263229, 1.6780178753252646}), 1e-6));
		EXPECT_TRUE(approxEqual(sol.u(1), (DynamicVector<double> {-2.5863408379530139}), 1e-6));
		//EXPECT_TRUE(approxEqual(sol.pi(e1), (DynamicVector<double> {37.522042463325405, 52.967530263604466}), 1e-6));

		EXPECT_TRUE(approxEqual(sol.x(2), (DynamicVector<double> {2.2238563940113898, 1.0916770373722506}), 1e-6));
	}


	TYPED_TEST_P(DynamicRiccatiTest, testSingleStageNoState)
	{
		size_t const N = 1, nx = 0, nu = 2;

        OcpTree const g(N + 1);
        DynamicOcpSize const sz {g, nx, nu, 0, 0, 0, true};
        
        typename TestFixture::Qp qp {sz};
		typename TestFixture::Solution sol {sz};
        typename TestFixture::RiccatiImpl riccati {sz};

        qp.R(vertex(0, g), blaze::DynamicMatrix<double> {
            {0.444923,    0.598522}, 
            {0.598522,     1.11451}
        });

        qp.r(vertex(0, g), blaze::DynamicVector<double> {
            0.131991, 0.585785
        });

        riccati(qp, sol);

        auto const v = vertex(0, g);
		EXPECT_TRUE(approxEqual(sol.u(v), blaze::DynamicVector<double> {1.4785, -1.3196}, 1e-4));
	}


	TYPED_TEST_P(DynamicRiccatiTest, testKktValue)
	{
		size_t constexpr NC = 0;
		size_t constexpr NCT = 0;

		OcpTree const g {
			3,
			2, 2, 2, 
			1, 1, 1, 1, 1, 1,
			0, 0, 0, 0, 0, 0
		};

		for (size_t nx = 1; nx < 10; ++nx)
		{
			for (size_t nu = 1; nu < 10; ++nu)
			{
				DynamicOcpSize const size {g, nx, nu, NC, 0, NCT, false};
				
				typename TestFixture::Qp qp {size};
				typename TestFixture::Solution sol {size};
				randomize(qp);
				
				typename TestFixture::RiccatiImpl riccati {size};
				riccati(qp, sol);

				// Set Lagrange multipliers for inequalities to 0
				for (auto u : vertices(g))
				{
					reset(sol.lam_lx(u));
					reset(sol.lam_ux(u));
					reset(sol.lam_ld(u));
					reset(sol.lam_ud(u));
				}

				for (auto u : g.branchVertices())
				{
					reset(sol.lam_lu(u));
					reset(sol.lam_uu(u));
				}

				EXPECT_TRUE(isSolution(sol, qp, 1e-10)) << "with nx=" << nx << ", nu=" << nu;
			}
		}
	}


    REGISTER_TYPED_TEST_SUITE_P(DynamicRiccatiTest
        , test_riccati0
		, test_riccati1
		, testSingleStageNoState
		, testKktValue
    );
}