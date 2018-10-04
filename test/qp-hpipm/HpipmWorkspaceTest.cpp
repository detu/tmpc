#include "TreeQpWorkspaceTest.hpp"
#include "TreeQpWorkspaceSolveTest.hpp"
#include "SoftConstraintsTest.hpp"

#include <tmpc/qp/HpipmWorkspace.hpp>
//#include <tmpc/EigenKernel.hpp>
#include <tmpc/BlazeKernel.hpp>


namespace tmpc :: testing
{
	INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, TreeQpWorkspaceTest, HpipmWorkspace<BlazeKernel<double>>);
	INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, TreeQpWorkspaceSolveTest, HpipmWorkspace<BlazeKernel<double>>);
	//INSTANTIATE_TYPED_TEST_CASE_P(Hpipm_Blaze_double, SoftConstraintsTest, HpipmWorkspace<BlazeKernel<double>>);


	TEST(HpipmWorkspaceTest, test_solveUnconstrained)
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
		using Workspace = HpipmWorkspace<BlazeKernel<double>>;
		using StageHessianMatrix = StaticMatrix<double, NZ, NZ>;

		OcpGraph const g = ocpGraphLinear(NT + 1);
		Workspace ws {g, ocpSizeNominalMpc(NT, NX, NU, NC, 0, NCT, false)};
		auto const e0 = out_edges(0, g).front();
		auto const e1 = out_edges(1, g).front();
		
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

		ws.solveUnconstrained();

		EXPECT_PRED2(MatrixApproxEquality(1e-6), get(ws.x(), 0), (DynamicVector<double> {4.2376727217537882, -3.2166970575479454}));
		EXPECT_PRED2(MatrixApproxEquality(1e-6), get(ws.u(), 0), (DynamicVector<double> {-0.34889319983994238}));

		EXPECT_PRED2(MatrixApproxEquality(1e-6), get(ws.x(), 1), (DynamicVector<double> {1.8465290642858716, -1.5655902573878877}));
		EXPECT_PRED2(MatrixApproxEquality(1e-6), get(ws.u(), 1), (DynamicVector<double> {-0.20287931724419153}));

		EXPECT_PRED2(MatrixApproxEquality(1e-6), get(ws.x(), 2), (DynamicVector<double> {1.1794991482758881, 0.23153042536792068}));
	}
}