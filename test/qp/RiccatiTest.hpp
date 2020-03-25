#pragma once

#include <tmpc/qp/MpipmWorkspace.hpp>
#include <tmpc/qp/OcpQp.hpp>
#include <tmpc/qp/KktResidual.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/Testing.hpp>
#include <tmpc/testing/qp/IsSolution.hpp>

#include <blaze/Math.h>


namespace tmpc :: testing
{
    template <typename RiccatiImpl_>
    class RiccatiTest 
    :   public Test
    {
    public:
        using Workspace = MpipmWorkspace<double>;
		using RiccatiImpl = RiccatiImpl_;

    protected:
    };


    TYPED_TEST_SUITE_P(RiccatiTest);


    TYPED_TEST_P(RiccatiTest, test_riccati0)
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
		using Workspace = typename TestFixture::Workspace;
		using StageHessianMatrix = StaticMatrix<double, NZ, NZ>;

		OcpGraph const g = ocpGraphLinear(NT + 1);
		auto const sz = ocpSizeNominalMpc(NT, NX, NU, NC, 0, NCT, false);
		Workspace ws {g, sz};

		auto const e0 = graph::out_edges(0, g).front();
		auto const e1 = graph::out_edges(1, g).front();
		
		put(ws.lx(), 0, DynamicVector<double>(NX, -1.));
		put(ws.lu(), 0, DynamicVector<double>(NU, -1.));
		put(ws.ux(), 0, DynamicVector<double>(NX, 1.));
		put(ws.uu(), 0, DynamicVector<double>(NU, 1.));

		put(ws.lx(), 1, DynamicVector<double>(NX, -1.));
		put(ws.lu(), 1, DynamicVector<double>(NU, -1.));
		put(ws.ux(), 1, DynamicVector<double>(NX, 1.));
		put(ws.uu(), 1, DynamicVector<double>(NU, 1.));

		put(ws.lx(), 2, DynamicVector<double>(NX, -1.));
		put(ws.ux(), 2, DynamicVector<double>(NX, 1.));
		
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
		put(ws.Q(), 0, Q0);	
		put(ws.R(), 0, R0);	
		put(ws.S(), 0, S0);	
		put(ws.q(), 0, q0);	
		put(ws.r(), 0, r0);

		put(ws.Q(), 1, Q1);	
		put(ws.R(), 1, R1);	
		put(ws.S(), 1, S1);	
		put(ws.q(), 1, q1);	
		put(ws.r(), 1, r1);

		put(ws.Q(), 2, Q2);	
		put(ws.q(), 2, q2);

		put(ws.A(), e0, A0);	
		put(ws.B(), e0, B0);	
		put(ws.b(), e0, a0);

		put(ws.A(), e1, A1);	
		put(ws.B(), e1, B1);	
		put(ws.b(), e1, a1);

		typename TestFixture::RiccatiImpl riccati(g, sz);
		riccati(ws, ws);

		EXPECT_TRUE(approxEqual(get(ws.x(), 0), (DynamicVector<double> {4.2376727217537882, -3.2166970575479454}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 0), (DynamicVector<double> {-0.34889319983994238}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.pi(), e0), (DynamicVector<double> {-1.621313883169287, 11.132830578594472}), 1e-6));

		EXPECT_TRUE(approxEqual(get(ws.x(), 1), (DynamicVector<double> {1.8465290642858716, -1.5655902573878877}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 1), (DynamicVector<double> {-0.20287931724419153}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.pi(), e1), (DynamicVector<double> {15.036417437909771, 21.143596583220845}), 1e-6));

		EXPECT_TRUE(approxEqual(get(ws.x(), 2), (DynamicVector<double> {1.1794991482758881, 0.23153042536792068}), 1e-6));
	}


	TYPED_TEST_P(RiccatiTest, test_riccati1)
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
		using Workspace = typename TestFixture::Workspace;
		using StageHessianMatrix = StaticMatrix<double, NZ, NZ>;

		OcpGraph const g = ocpGraphLinear(NT + 1);
		auto const sz = ocpSizeNominalMpc(NT, NX, NU, NC, 0, NCT, true);
		Workspace ws {g, sz};

		auto const e0 = graph::out_edges(0, g).front();
		auto const e1 = graph::out_edges(1, g).front();
		
		put(ws.lx(), 0, DynamicVector<double>(0, -1.));
		put(ws.lu(), 0, DynamicVector<double>(NU, -1.));
		put(ws.ux(), 0, DynamicVector<double>(0, 1.));
		put(ws.uu(), 0, DynamicVector<double>(NU, 1.));

		put(ws.lx(), 1, DynamicVector<double>(NX, -1.));
		put(ws.lu(), 1, DynamicVector<double>(NU, -1.));
		put(ws.ux(), 1, DynamicVector<double>(NX, 1.));
		put(ws.uu(), 1, DynamicVector<double>(NU, 1.));

		put(ws.lx(), 2, DynamicVector<double>(NX, -1.));
		put(ws.ux(), 2, DynamicVector<double>(NX, 1.));
		
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
		put(ws.Q(), 0, Q0);	
		put(ws.R(), 0, R0);	
		put(ws.S(), 0, S0);	
		put(ws.q(), 0, q0);	
		put(ws.r(), 0, r0);

		put(ws.Q(), 1, Q1);	
		put(ws.R(), 1, R1);	
		put(ws.S(), 1, S1);	
		put(ws.q(), 1, q1);	
		put(ws.r(), 1, r1);

		put(ws.Q(), 2, Q2);	
		put(ws.q(), 2, q2);

		put(ws.A(), e0, A0);	
		put(ws.B(), e0, B0);	
		put(ws.b(), e0, a0);

		put(ws.A(), e1, A1);	
		put(ws.B(), e1, B1);	
		put(ws.b(), e1, a1);

		typename TestFixture::RiccatiImpl riccati(g, sz);
		riccati(ws, ws);

		EXPECT_TRUE(approxEqual(get(ws.x(), 0), (DynamicVector<double> {}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 0), (DynamicVector<double> {-0.32198212467473536}), 1e-6));
		//EXPECT_TRUE(approxEqual(get(ws.pi(), e0), (DynamicVector<double> {-8.1496398536788472, 34.341139646264558}), 1e-6));

		EXPECT_TRUE(approxEqual(get(ws.x(), 1), (DynamicVector<double> {0.83900893766263229, 1.6780178753252646}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 1), (DynamicVector<double> {-2.5863408379530139}), 1e-6));
		//EXPECT_TRUE(approxEqual(get(ws.pi(), e1), (DynamicVector<double> {37.522042463325405, 52.967530263604466}), 1e-6));

		EXPECT_TRUE(approxEqual(get(ws.x(), 2), (DynamicVector<double> {2.2238563940113898, 1.0916770373722506}), 1e-6));
	}


	TYPED_TEST_P(RiccatiTest, testSingleStageNoState)
	{
		size_t const N = 1, nx = 0, nu = 2;

        OcpGraph const g = ocpGraphLinear(N + 1);
        auto const sz = ocpSizeNominalMpc(N, nx, nu, 0, 0, 0, true);
        
        MpipmWorkspace<double> ws_mpipm(g, sz);
        typename TestFixture::RiccatiImpl riccati(g, sz);

        put(ws_mpipm.R(), vertex(0, g), blaze::DynamicMatrix<double> {
            {0.444923,    0.598522}, 
            {0.598522,     1.11451}
        });

        put(ws_mpipm.r(), vertex(0, g), blaze::DynamicVector<double> {
            0.131991, 0.585785
        });

        riccati(ws_mpipm, ws_mpipm);

        auto const v = vertex(0, g);
		EXPECT_TRUE(approxEqual(get(ws_mpipm.u(), v), blaze::DynamicVector<double> {1.4785, -1.3196}, 1e-4));
	}


	TYPED_TEST_P(RiccatiTest, testKktResidual)
	{
		using Real = double;
		using Workspace = typename TestFixture::Workspace;

		size_t constexpr NX = 2;
		size_t constexpr NU = 1;
		size_t constexpr NC = 0;
		size_t constexpr NCT = 0;
		size_t constexpr NT = 2;

		OcpGraph const g = ocpGraphLinear(NT + 1);
		auto const size_map = ocpSizeNominalMpc(NT, NX, NU, NC, 0, NCT, false);
		
		MpipmWorkspace<Real> ws(g, size_map);
        randomizeQp(ws);
		
		typename TestFixture::RiccatiImpl riccati(g, size_map);
		riccati(ws, ws);

		EXPECT_TRUE(isSolution(ws, ws, 1e-14));
	}


    REGISTER_TYPED_TEST_SUITE_P(RiccatiTest
        , test_riccati0
		, test_riccati1
		, testSingleStageNoState
		, testKktResidual
    );
}