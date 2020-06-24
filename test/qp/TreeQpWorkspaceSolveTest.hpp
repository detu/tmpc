#pragma once

#include <tmpc/Matrix.hpp>
#include <tmpc/Math.hpp>
#include <tmpc/ocp/OcpGraph.hpp>
#include <tmpc/ocp/OcpSizeProperties.hpp>
#include <tmpc/Traits.hpp>

#include <tmpc/Testing.hpp>

#include <iostream>
#include <fstream>


namespace tmpc :: testing
{
	template <typename WS>
	class TreeQpWorkspaceSolveTest 
	: 	public Test
	{
	protected:
		using Workspace = WS;
		using Real = typename RealOf<WS>::type;
		using Vector = blaze::DynamicVector<Real, blaze::columnVector>;
		using Matrix = blaze::DynamicMatrix<Real>;

		static Workspace problem_0()
		{
			unsigned const NX = 2;
			unsigned const NU = 1;
			unsigned const NZ = NX + NU;
			unsigned const NC = 0;
			unsigned const NCT = 0;
			unsigned const NT = 2;

			using StageHessianMatrix = blaze::StaticMatrix<Real, NZ, NZ>;

			OcpGraph const g = ocpGraphLinear(NT + 1);
			Workspace ws {g, ocpSizeNominalMpc(NT, NX, NU, NC, 0, NCT, false)};
			auto const e0 = graph::out_edges(0, g).front();
			auto const e1 = graph::out_edges(1, g).front();

			/*
			auto qp = ws.problem();
			qp[0].lbx(-1.);	qp[0].lbu(-1.);	qp[0].ubx(1.);	qp[0].ubu(1.);
			qp[1].lbx(-1.);	qp[1].lbu(-1.);	qp[1].ubx(1.);	qp[1].ubu(1.);
			qp[2].lbx(-1.);					qp[2].ubx(1.);
			*/
			
			put(ws.lx(), 0, Vector(NX, -1.));
			put(ws.lu(), 0, Vector(NU, -1.));
			put(ws.ux(), 0, Vector(NX, 1.));
			put(ws.uu(), 0, Vector(NU, 1.));

			put(ws.lx(), 1, Vector(NX, -1.));
			put(ws.lu(), 1, Vector(NU, -1.));
			put(ws.ux(), 1, Vector(NX, 1.));
			put(ws.uu(), 1, Vector(NU, 1.));

			put(ws.lx(), 2, Vector(NX, -1.));
			put(ws.ux(), 2, Vector(NX, 1.));
			
			// Stage 0
			StageHessianMatrix H0 {
				{67,   78,   90},
				{78,   94,  108},
				{90,  108,  127}
			};

			blaze::StaticVector<Real, NX> const q0 {0., 0.};
			blaze::StaticVector<Real, NU> const r0 {0.};

			const Matrix Q0 = submatrix(H0, 0, 0, NX, NX);
			const Matrix R0 = submatrix(H0, NX, NX, NU, NU);
			const Matrix S0T = submatrix(H0, 0, NX, NX, NU);
			const Matrix S0 = submatrix(H0, NX, 0, NU, NX);

			Matrix const A0 {{1., 1.}, {0., 1.}};
			Matrix const B0 {{0.5}, {1.0}};
			Vector a0 {1., 2.};

			// Stage 1
			StageHessianMatrix H1 {
				{67,   78,   90},
				{78,   94,  108},
				{90,  108,  127}
			};

			blaze::StaticVector<Real, NX> const q1 {0., 0.};
			blaze::StaticVector<Real, NU> const r1 {0.};

			const Matrix Q1 = submatrix(H1, 0, 0, NX, NX);
			const Matrix R1 = submatrix(H1, NX, NX, NU, NU);
			const Matrix S1T = submatrix(H1, 0, NX, NX, NU);
			const Matrix S1 = submatrix(H1, NX, 0, NU, NX);

			Matrix const A1 {{1., 1.}, {0., 1.}};
			Matrix const B1 {{0.5}, {1.0}};
			Vector const a1 {1., 2.};

			// Stage 2
			Matrix H2 {{1., 2.}, {3., 4.}};
			H2 = trans(H2) * H2;	// Make positive definite.

			blaze::StaticVector<Real, NX> const q2 {0., 0.};

			const Matrix Q2 = submatrix(H2, 0, 0, NX, NX);

			// Setup QP
			/*
			qp[0].Q(Q0);	qp[0].R(R0);	qp[0].S(S0);	qp[0].q(q0);	qp[0].r(r0);
			qp[1].Q(Q1);	qp[1].R(R1);	qp[1].S(S1);	qp[1].q(q1);	qp[1].r(r1);
			qp[2].Q(Q2);									qp[2].q(q2);

			qp[0].A(A0);	
			qp[0].B(B0);		
			qp[0].b(a0);
			
			qp[1].A(A1);	
			qp[1].B(B1);		
			qp[1].b(a1);
			*/
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

			return std::move(ws);
		}


		static Workspace problem_0_unconstrained()
		{
			unsigned const NX = 2;
			unsigned const NU = 1;
			unsigned const NZ = NX + NU;
			unsigned const NC = 0;
			unsigned const NCT = 0;
			unsigned const NT = 2;

			typedef blaze::StaticMatrix<Real, NZ, NZ> StageHessianMatrix;

			OcpGraph const g = ocpGraphLinear(NT + 1);
			Workspace ws {g, ocpSizeNominalMpc(NT, NX, NU, NC, 0, NCT, false)};
			auto const e0 = graph::out_edges(0, g).front();
			auto const e1 = graph::out_edges(1, g).front();

			/*
			auto qp = ws.problem();
			qp[0].lbx(-1.);	qp[0].lbu(-1.);	qp[0].ubx(1.);	qp[0].ubu(1.);
			qp[1].lbx(-1.);	qp[1].lbu(-1.);	qp[1].ubx(1.);	qp[1].ubu(1.);
			qp[2].lbx(-1.);					qp[2].ubx(1.);
			*/
			
			put(ws.lx(), 0, Vector(NX, -inf<Real>()));
			put(ws.lu(), 0, Vector(NU, -inf<Real>()));
			put(ws.ux(), 0, Vector(NX, inf<Real>()));
			put(ws.uu(), 0, Vector(NU, inf<Real>()));

			put(ws.lx(), 1, Vector(NX, -inf<Real>()));
			put(ws.lu(), 1, Vector(NU, -inf<Real>()));
			put(ws.ux(), 1, Vector(NX, inf<Real>()));
			put(ws.uu(), 1, Vector(NU, inf<Real>()));

			put(ws.lx(), 2, Vector(NX, -inf<Real>()));
			put(ws.ux(), 2, Vector(NX, inf<Real>()));
			
			// Stage 0
			StageHessianMatrix H0 {
				{67,   78,   90},
				{78,   94,  108},
				{90,  108,  127}
			};

			blaze::StaticVector<Real, NX> const q0 {0., 0.};
			blaze::StaticVector<Real, NU> const r0 {0.};

			const Matrix Q0 = submatrix(H0, 0, 0, NX, NX);
			const Matrix R0 = submatrix(H0, NX, NX, NU, NU);
			const Matrix S0T = submatrix(H0, 0, NX, NX, NU);
			const Matrix S0 = submatrix(H0, NX, 0, NU, NX);

			Matrix const A0 {{1., 1.}, {0., 1.}};
			Matrix const B0 {{0.5}, {1.0}};
			Vector a0 {1., 2.};

			// Stage 1
			StageHessianMatrix H1 {
				{67,   78,   90},
				{78,   94,  108},
				{90,  108,  127}
			};

			blaze::StaticVector<Real, NX> const q1 {0., 0.};
			blaze::StaticVector<Real, NU> const r1 {0.};

			const Matrix Q1 = submatrix(H1, 0, 0, NX, NX);
			const Matrix R1 = submatrix(H1, NX, NX, NU, NU);
			const Matrix S1T = submatrix(H1, 0, NX, NX, NU);
			const Matrix S1 = submatrix(H1, NX, 0, NU, NX);

			Matrix const A1 {{1., 1.}, {0., 1.}};
			Matrix const B1 {{0.5}, {1.0}};
			Vector const a1 {1., 2.};

			// Stage 2
			Matrix H2 {{1., 2.}, {3., 4.}};
			H2 = trans(H2) * H2;	// Make positive definite.

			blaze::StaticVector<Real, NX> const q2 {0., 0.};

			const Matrix Q2 = submatrix(H2, 0, 0, NX, NX);

			// Setup QP
			/*
			qp[0].Q(Q0);	qp[0].R(R0);	qp[0].S(S0);	qp[0].q(q0);	qp[0].r(r0);
			qp[1].Q(Q1);	qp[1].R(R1);	qp[1].S(S1);	qp[1].q(q1);	qp[1].r(r1);
			qp[2].Q(Q2);									qp[2].q(q2);

			qp[0].A(A0);	
			qp[0].B(B0);		
			qp[0].b(a0);
			
			qp[1].A(A1);	
			qp[1].B(B1);		
			qp[1].b(a1);
			*/
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

			return std::move(ws);
		}


		/*
		static Workspace problem_1()
		{
			Workspace ws = problem_0();
			auto qp = ws.problem();

			blaze::StaticVector<Real, 2> x0 {1., 0.};
			qp[0].lbx(x0);	qp[0].ubx(x0);

			qp[0].b(0.);
			qp[1].b(0.);

			return std::move(ws);
		}
		*/
	};

	TYPED_TEST_SUITE_P(TreeQpWorkspaceSolveTest);

	/// \brief Check if QPSolver move constructor works and the solver works after move constructor.
	/*
	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testMoveConstructor)
	{
		auto ws = TestFixture::problem_0();
		auto ws1 = std::move(ws);

		ws1.solve();
		auto const sol = ws1.solution();

		EXPECT_PRED2(ApproxEqual(1e-6), sol[0].x(), (typename TestFixture::Vector {1., -1.}));
		EXPECT_PRED2(ApproxEqual(1e-6), sol[0].u(), (typename TestFixture::Vector {-1.}));

		EXPECT_PRED2(ApproxEqual(1e-6), sol[1].x(), (typename TestFixture::Vector {0.5, 0.}));
		EXPECT_PRED2(ApproxEqual(1e-6), sol[1].u(), (typename TestFixture::Vector {-1.}));

		EXPECT_PRED2(ApproxEqual(1e-6), sol[2].x(), (typename TestFixture::Vector {1., 1.}));
	}
	*/

	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testSolve0)
	{
		auto ws = TestFixture::problem_0();
		ws.solve();

		EXPECT_TRUE(approxEqual(get(ws.x(), 0), (typename TestFixture::Vector {1., -1.}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 0), (typename TestFixture::Vector {-1.}), 1e-6));

		EXPECT_TRUE(approxEqual(get(ws.x(), 1), (typename TestFixture::Vector {0.5, 0.}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 1), (typename TestFixture::Vector {-1.}), 1e-6));

		EXPECT_TRUE(approxEqual(get(ws.x(), 2), (typename TestFixture::Vector {1., 1.}), 1e-6));
	}


	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testSolve0Unconstrained)
	{
		auto ws = TestFixture::problem_0_unconstrained();
		ws.solve();

		EXPECT_TRUE(approxEqual(get(ws.x(), 0), (typename TestFixture::Vector {4.2376727217537882, -3.2166970575479454}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 0), (typename TestFixture::Vector {-0.34889319983994238}), 1e-6));

		EXPECT_TRUE(approxEqual(get(ws.x(), 1), (typename TestFixture::Vector {1.8465290642858716, -1.5655902573878877}), 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 1), (typename TestFixture::Vector {-0.20287931724419153}), 1e-6));

		EXPECT_TRUE(approxEqual(get(ws.x(), 2), (typename TestFixture::Vector {1.1794991482758881, 0.23153042536792068}), 1e-6));
	}


	// Build a small tree
    //
    //     1 - 3
    //   /
    // 0 - 2 - 4
	//
	// /* Dimensions */
	//
	// int Nn = 5;
	// int ns[5] = { 2, 1, 1, 0, 0, };
	// int nx[5] = { 2, 2, 2, 1, 4, };
	// int nu[5] = { 2, 2, 2, 0, 0, };
	//
	// /* Data */
	//
	// double A[18] = { 1.844336677576532e-01, 2.120308425323207e-01, 7.734680811267680e-02, 9.138004107795679e-01, 3.353113307052461e-01, 2.992250233331066e-01, 4.525925415693240e-01, 4.226456532204624e-01, 6.659872164111106e-01, 8.943893753542428e-01, 3.656301804845286e-02, 8.092038512937934e-01, 7.486188717761971e-01, 1.201870179870806e-01, 5.250451647626088e-01, 3.258336287632492e-01, 5.464494399030685e-01, 3.988807523831990e-01, };
	// double B[18] = { 7.067152176969306e-01, 5.577889667548762e-01, 3.134289899365913e-01, 1.662035629021507e-01, 3.596063179722356e-01, 5.583191998692971e-01, 7.425453657019391e-01, 4.243347836256907e-01, 5.165582083512704e-01, 7.027023069504753e-01, 4.150933866130466e-01, 1.807377602547944e-01, 2.553867404880508e-01, 2.053577465818457e-02, 9.236756126204072e-01, 6.536998890082529e-01, 9.326135720485641e-01, 1.635123685275256e-01, };
	// double b[9] = { 6.224972592798952e-01, 9.879347349524954e-01, 4.293557885762050e-01, 1.248727587198128e-01, 1.535903766194002e-01, 9.210972558921975e-01, 7.946578853887531e-01, 5.773941967066487e-01, 4.400355957602536e-01, };
	// double Q[29] = { 5.319913229393722e+00, 5.393412800510404e-01, 5.393412800510404e-01, 1.811101907422479e+00, 7.011392692676973e+00, 3.272957846025784e-01, 3.272957846025784e-01, 4.097878096538557e+00, 9.593793256757216e+00, 3.038529240149768e-01, 3.038529240149768e-01, 1.001099886181528e+01, 6.362297882301013e+00, 6.614227625326206e+00, 7.596379523220123e-01, 3.238588878651922e-01, 4.393634301462749e-01, 7.596379523220123e-01, 1.018014633913789e+01, 5.529872984950073e-01, 7.155671716864053e-01, 3.238588878651922e-01, 5.529872984950073e-01, 5.255781653063426e+00, 5.848458769998477e-01, 4.393634301462749e-01, 7.155671716864053e-01, 5.848458769998477e-01, 1.105824815366516e+00, };
	// double R[12] = { 9.158208703923456e-01, 1.304015848214388e-01, 1.304015848214388e-01, 4.350966896940532e-01, 1.745525610282895e+00, 5.114279662110315e-01, 5.114279662110315e-01, 4.766377491779941e-01, 1.939182463811595e+00, 5.021881705592814e-01, 5.021881705592814e-01, 2.246704066105748e+00, };
	// double S[12] = { 1.962489222569553e-01, 3.174797751494354e-01, 3.164289991462910e-01, 2.175633094228206e-01, 7.581124313274186e-01, 8.711111219153892e-01, 3.507767448858926e-01, 6.855357087475372e-01, 1.059204167327653e-01, 6.815604304703157e-01, 4.632605785937192e-01, 2.121632052549344e-01, };
	// double q[11] = { 2.510418460157361e-01, 8.929224052859770e-01, 2.941486337678496e-01, 5.306293038568856e-01, 9.851873768810837e-02, 8.235744739278386e-01, 6.797338982104669e-01, 8.667498969993187e-01, 6.311887342690112e-01, 3.550736518788490e-01, 9.970032716066477e-01, };
	// double r[6] = { 7.032232245562910e-01, 5.557379427193866e-01, 8.324233862851839e-01, 5.974901918725793e-01, 1.750097373820796e-01, 1.635699097849932e-01, };
	//
	// /* Optimal solution */
	//
	// double xopt[11] = { 6.592835368068550e-02, -4.063901103087715e-01, 2.342956460060519e-01, 3.561800465271370e-01, -4.093033514742150e-02, -3.552137307936658e-01, -2.087478975092527e-01, 2.939583174994020e-02, 1.623760005176457e-01, -3.367157300712413e-01, 1.756957363329313e-01, };
	// double uopt[6] = { -4.301938159361802e-01, -2.070756990391132e-01, -1.394097954236166e-01, -1.088549940479685e+00, -1.291731382101928e-01, -7.037998096601095e-01, };
	//
	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testSolveDimitrisRandom)
	{
		using Vec = typename TestFixture::Vector;
		using Mat = typename TestFixture::Matrix;
		using Real = typename TestFixture::Real;

		size_t const N = 5;
		std::array<size_t, N> const n_kids = {2, 1, 1, 0, 0};
		std::array<OcpSize, N> const sz = {
			OcpSize {2, 2},
			OcpSize {2, 2},
			OcpSize {2, 2},
			OcpSize {1, 0},
			OcpSize {4, 0}
		};

		OcpGraph const g = ocpGraphFromOutDegree(n_kids.begin());
		typename TestFixture::Workspace ws {g, sz.begin()};

		for (auto v : boost::make_iterator_range(vertices(g)))
		{
			put(ws.lx(), v, Vec(sz[v].nx(), -inf<Real>()));
			put(ws.lu(), v, Vec(sz[v].nu(), -inf<Real>()));
			put(ws.ux(), v, Vec(sz[v].nx(), inf<Real>()));
			put(ws.uu(), v, Vec(sz[v].nu(), inf<Real>()));
		}

		auto const v0 = vertex(0, g);
		auto const v1 = vertex(1, g);
		auto const v2 = vertex(2, g);
		auto const v3 = vertex(3, g);
		auto const v4 = vertex(4, g);

		put(ws.Q(), v0, Mat {
			{5.319913229393722e+00, 5.393412800510404e-01}, 
			{5.393412800510404e-01, 1.811101907422479e+00}
		});

		put(ws.Q(), v1, Mat {
			{7.011392692676973e+00, 3.272957846025784e-01}, 
			{3.272957846025784e-01, 4.097878096538557e+00}
		});

		put(ws.Q(), v2, Mat {
			{9.593793256757216e+00, 3.038529240149768e-01}, 
			{3.038529240149768e-01, 1.001099886181528e+01}
		});

		put(ws.Q(), v3, Mat {
			{6.362297882301013e+00}
		});

		put(ws.Q(), v4, Mat {
			{6.614227625326206e+00, 7.596379523220123e-01, 3.238588878651922e-01, 4.393634301462749e-01},
			{7.596379523220123e-01, 1.018014633913789e+01, 5.529872984950073e-01, 7.155671716864053e-01}, 
			{3.238588878651922e-01, 5.529872984950073e-01, 5.255781653063426e+00, 5.848458769998477e-01}, 
			{4.393634301462749e-01, 7.155671716864053e-01, 5.848458769998477e-01, 1.105824815366516e+00}
		});

		put(ws.R(), v0, Mat {
			{9.158208703923456e-01, 1.304015848214388e-01}, 
			{1.304015848214388e-01, 4.350966896940532e-01}
		});

		put(ws.R(), v1, Mat {
			{1.745525610282895e+00, 5.114279662110315e-01}, 
			{5.114279662110315e-01, 4.766377491779941e-01}
		});

		put(ws.R(), v2, Mat {
			{1.939182463811595e+00, 5.021881705592814e-01}, 
			{5.021881705592814e-01, 2.246704066105748e+00}
		});

		put(ws.S(), v0, trans(Mat {
			{1.962489222569553e-01, 3.174797751494354e-01}, 
			{3.164289991462910e-01, 2.175633094228206e-01}
		}));

		put(ws.S(), v1, trans(Mat {
			{7.581124313274186e-01, 8.711111219153892e-01}, 
			{3.507767448858926e-01, 6.855357087475372e-01}
		}));

		put(ws.S(), v2, trans(Mat {
			{1.059204167327653e-01, 6.815604304703157e-01}, 
			{4.632605785937192e-01, 2.121632052549344e-01}
		}));

		put(ws.q(), v0, Vec {2.510418460157361e-01, 8.929224052859770e-01});
		put(ws.q(), v1, Vec {2.941486337678496e-01, 5.306293038568856e-01});
		put(ws.q(), v2, Vec {9.851873768810837e-02, 8.235744739278386e-01});
		put(ws.q(), v3, Vec {6.797338982104669e-01});
		put(ws.q(), v4, Vec {8.667498969993187e-01, 6.311887342690112e-01, 3.550736518788490e-01, 9.970032716066477e-01});

		put(ws.r(), v0, Vec {7.032232245562910e-01, 5.557379427193866e-01});
		put(ws.r(), v1, Vec {8.324233862851839e-01, 5.974901918725793e-01});
		put(ws.r(), v2, Vec {1.750097373820796e-01, 1.635699097849932e-01});

		auto const e1 = edge(v0, v1, g).first;
		auto const e2 = edge(v0, v2, g).first;
		auto const e3 = edge(v1, v3, g).first;
		auto const e4 = edge(v2, v4, g).first;

		put(ws.A(), e1, trans(Mat {
			{1.844336677576532e-01, 2.120308425323207e-01}, 
			{7.734680811267680e-02, 9.138004107795679e-01}
		}));

		put(ws.A(), e2, trans(Mat {
			{3.353113307052461e-01, 2.992250233331066e-01}, 
			{4.525925415693240e-01, 4.226456532204624e-01}
		}));

		put(ws.A(), e3, trans(Mat {
			{6.659872164111106e-01}, 
			{8.943893753542428e-01}
		}));

		put(ws.A(), e4, trans(Mat {
			{3.656301804845286e-02, 8.092038512937934e-01, 7.486188717761971e-01, 1.201870179870806e-01}, 
			{5.250451647626088e-01, 3.258336287632492e-01, 5.464494399030685e-01, 3.988807523831990e-01}
		}));

		put(ws.B(), e1, trans(Mat {
			{7.067152176969306e-01, 5.577889667548762e-01}, 
			{3.134289899365913e-01, 1.662035629021507e-01}
		}));

		put(ws.B(), e2, trans(Mat {
			{3.596063179722356e-01, 5.583191998692971e-01}, 
			{7.425453657019391e-01, 4.243347836256907e-01}
		}));

		put(ws.B(), e3, trans(Mat {
			{5.165582083512704e-01}, 
			{7.027023069504753e-01}
		}));

		put(ws.B(), e4, trans(Mat {
			{4.150933866130466e-01, 1.807377602547944e-01, 2.553867404880508e-01, 2.053577465818457e-02}, 
			{9.236756126204072e-01, 6.536998890082529e-01, 9.326135720485641e-01, 1.635123685275256e-01}
		}));

		put(ws.b(), e1, Vec {6.224972592798952e-01, 9.879347349524954e-01});
		put(ws.b(), e2, Vec {4.293557885762050e-01, 1.248727587198128e-01});
		put(ws.b(), e3, Vec {1.535903766194002e-01});
		put(ws.b(), e4, Vec {9.210972558921975e-01, 7.946578853887531e-01, 5.773941967066487e-01, 4.400355957602536e-01});

		ws.solve();

		EXPECT_TRUE(approxEqual(get(ws.x(), 0), Vec {6.592835368068550e-02, -4.063901103087715e-01}, 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.x(), 1), Vec {2.342956460060519e-01, 3.561800465271370e-01}, 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.x(), 2), Vec {-4.093033514742150e-02, -3.552137307936658e-01}, 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.x(), 3), Vec {-2.087478975092527e-01}, 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.x(), 4), Vec {2.939583174994020e-02, 1.623760005176457e-01, -3.367157300712413e-01, 1.756957363329313e-01}, 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 0), Vec {-4.301938159361802e-01, -2.070756990391132e-01}, 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 1), Vec {-1.394097954236166e-01, -1.088549940479685e+00}, 1e-6));
		EXPECT_TRUE(approxEqual(get(ws.u(), 2), Vec {-1.291731382101928e-01, -7.037998096601095e-01}, 1e-6));
	}

#if 0
	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testSolve1)
	{
		auto ws = TestFixture::problem_1();

		ws.solve();
		auto const sol = ws.solution();

		EXPECT_PRED2(ApproxEqual(1e-6), sol[0].x(), (typename TestFixture::Vector {1., 0.}));
		EXPECT_PRED2(ApproxEqual(1e-6), sol[0].u(), (typename TestFixture::Vector {-0.690877362606266}));

		EXPECT_PRED2(ApproxEqual(1e-6), sol[1].x(), (typename TestFixture::Vector {0.654561318696867, -0.690877362606266}));
		EXPECT_PRED2(ApproxEqual(1e-6), sol[1].u(), (typename TestFixture::Vector {0.215679569867116}));

		EXPECT_PRED2(ApproxEqual(1e-6), sol[2].x(), (typename TestFixture::Vector {0.0715237410241597, -0.475197792739149}));
	}

	TYPED_TEST_P(TreeQpWorkspaceSolveTest, DISABLED_testSolve1stage1d)
	{
		typename TestFixture::Workspace ws { std::array<OcpSize, 1> { OcpSize(1, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(2.);
		problem[0].q(-1.);
		problem[0].lbx(-100.);
		problem[0].ubx(100.);

		ws.solve();
		auto solution = ws.solution();

		EXPECT_PRED2(ApproxEqual(1e-6), solution[0].x(), (typename TestFixture::Vector {0.5}));
	}

	///
	/// minimize (1/2 * Q0 * x0^2 - q0 * x0) + (1/2 * Q1 * x1^2 - q1 * x1)
	/// s.t. x1 = A0 * x0 + b0
	///      -100 <= x0 <= 100
	///      -100 <= x1 <= 100
	///
	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testSolve2stage1d)
	{
		using Real = typename TestFixture::Real;

		Real const Q0 = 2., q0 = -1., Q1 = 1., q1 = 0., A0 = 0.2, b0 = 1.;

		typename TestFixture::Workspace ws {std::array<OcpSize, 2> { OcpSize(1, 0, 0), OcpSize(1, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(Q0);
		problem[0].q(q0);
		problem[0].A(A0);
		problem[0].b(b0);
		problem[0].lbx(-100.);
		problem[0].ubx(100.);

		problem[1].Q(Q1);
		problem[1].q(q1);
		problem[1].lbx(-100.);
		problem[1].ubx(100.);

		ws.solve();
		auto solution = ws.solution();

		Real const x0_opt = -(q0 + Q1 * A0 * b0 + q1 * A0) / (Q0 + Q1 * A0 * A0);

		EXPECT_PRED2(ApproxEqual(1e-6), solution[0].x(), (typename TestFixture::Vector {x0_opt}));
		EXPECT_PRED2(ApproxEqual(1e-6), solution[1].x(), (typename TestFixture::Vector {A0 * x0_opt + b0}));
	}

	/**
	 \brief This QP has 3 states and 2 inputs, so that the S matrix is 3x2.
	 
	 It checks if the element order of S (row-major or col-major) is correct 
	 for the solvers that depend on it (like hpmpc or hpipm).
	 */
	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testSolve2stage5d)
	{
		using Real = typename TestFixture::Real;

		typename TestFixture::Workspace ws {std::array<OcpSize, 2> { OcpSize(3, 2, 0), OcpSize(0, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(DynamicMatrix<typename TestFixture::Kernel>({
			{2.73449304081811,	1.88589647487061,	2.07850896302612},
			{1.88589647487061,	2.23398400982993,	2.04607087035366},
			{2.07850896302612,	2.04607087035366,	2.75909914758841}
		}));
		
		problem[0].R(DynamicMatrix<typename TestFixture::Kernel>({
			{2.58480403401811,	2.27683785836085},
			{2.27683785836085,	2.48531656389865}
		}));

		problem[0].S(DynamicMatrix<typename TestFixture::Kernel>({
			{1.94423942432824,	1.95671439063179},
			{2.31644542710892,	2.08748866537729},
			{2.46061457306194,	1.94728938270016}
		}));

		problem[0].q(typename TestFixture::Vector({
			0.276025076998578,
			0.679702676853675,
			0.655098003973841
		}));

		problem[0].r(typename TestFixture::Vector({
			0.162611735194631,
			0.118997681558377
		}));

		problem[0].A(StaticMatrix<typename TestFixture::Kernel, 0, 3>(0.));
		problem[0].B(StaticMatrix<typename TestFixture::Kernel, 0, 2>(0.));
		problem[0].b(StaticVector<typename TestFixture::Kernel, 0>(0.));
		problem[0].lbx(-10000.);
		problem[0].ubx(10000.);
		problem[0].lbu(-10000.);
		problem[0].ubu(10000.);

		ws.solve();
		auto solution = ws.solution();

		EXPECT_PRED2(ApproxEqual(1e-6), solution[0].x(), (typename TestFixture::Vector {
			146.566682434017,
			-427.218558989821,
			-345.347969700289
		}));

		EXPECT_PRED2(ApproxEqual(1e-6), solution[0].u(), (typename TestFixture::Vector {
			769.663140469139,
			-191.122524763114
		}));

		EXPECT_PRED2(ApproxEqual(1e-6), solution[1].x(), (typename TestFixture::Vector {}));
		EXPECT_PRED2(ApproxEqual(1e-6), solution[1].u(), (typename TestFixture::Vector {}));
	}

	/**
	 \brief Check that the case when all state and input bounds are infinite is correctly handled by a QP solver.
	 */
	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testInfiniteBoundsAll)
	{
		using Real = typename TestFixture::Real;

		typename TestFixture::Workspace ws {std::array<OcpSize, 2> { OcpSize(3, 2, 0), OcpSize(0, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(DynamicMatrix<typename TestFixture::Kernel>({
			{2.73449304081811,	1.88589647487061,	2.07850896302612},
			{1.88589647487061,	2.23398400982993,	2.04607087035366},
			{2.07850896302612,	2.04607087035366,	2.75909914758841}
		}));
		
		problem[0].R(DynamicMatrix<typename TestFixture::Kernel>({
			{2.58480403401811,	2.27683785836085},
			{2.27683785836085,	2.48531656389865}
		}));

		problem[0].S(DynamicMatrix<typename TestFixture::Kernel>({
			{1.94423942432824,	1.95671439063179},
			{2.31644542710892,	2.08748866537729},
			{2.46061457306194,	1.94728938270016}
		}));

		problem[0].q(typename TestFixture::Vector({
			0.276025076998578,
			0.679702676853675,
			0.655098003973841
		}));

		problem[0].r(typename TestFixture::Vector({
			0.162611735194631,
			0.118997681558377
		}));

		problem[0].A(StaticMatrix<typename TestFixture::Kernel, 0, 3>(0.));
		problem[0].B(StaticMatrix<typename TestFixture::Kernel, 0, 2>(0.));
		problem[0].b(StaticVector<typename TestFixture::Kernel, 0>(0.));
		problem[0].lbx(-inf<Real>());
		problem[0].ubx(inf<Real>());
		problem[0].lbu(-inf<Real>());
		problem[0].ubu(inf<Real>());

		ws.solve();
		auto solution = ws.solution();

		using Vector = typename TestFixture::Vector;

		EXPECT_PRED2(ApproxEqual(1e-6), solution[0].x(), (Vector {
			146.566682434017,
			-427.218558989821,
			-345.347969700289
		}));

		EXPECT_PRED2(ApproxEqual(1e-6), solution[0].u(), (Vector {
			769.663140469139,
			-191.122524763114
		}));

		EXPECT_PRED2(ApproxEqual(1e-6), solution[1].x(), (Vector {}));
		EXPECT_PRED2(ApproxEqual(1e-6), solution[1].u(), (Vector {}));
	}

	/**
	 \brief Check that the case when some state and input bounds are infinite is correctly handled by a QP solver.
	 */
	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testInfiniteBoundsSome)
	{
		using Real = typename TestFixture::Real;

		typename TestFixture::Workspace ws {std::array<OcpSize, 2> { OcpSize(3, 2, 0), OcpSize(0, 0, 0) } };

		auto problem = ws.problem();
		problem[0].Q(DynamicMatrix<typename TestFixture::Kernel>({
			{2.73449304081811,	1.88589647487061,	2.07850896302612},
			{1.88589647487061,	2.23398400982993,	2.04607087035366},
			{2.07850896302612,	2.04607087035366,	2.75909914758841}
		}));
		
		problem[0].R(DynamicMatrix<typename TestFixture::Kernel>({
			{2.58480403401811,	2.27683785836085},
			{2.27683785836085,	2.48531656389865}
		}));

		problem[0].S(DynamicMatrix<typename TestFixture::Kernel>({
			{1.94423942432824,	1.95671439063179},
			{2.31644542710892,	2.08748866537729},
			{2.46061457306194,	1.94728938270016}
		}));

		problem[0].q(typename TestFixture::Vector({
			0.276025076998578,
			0.679702676853675,
			0.655098003973841
		}));

		problem[0].r(typename TestFixture::Vector({
			0.162611735194631,
			0.118997681558377
		}));

		problem[0].A(StaticMatrix<typename TestFixture::Kernel, 0, 3>(0.));
		problem[0].B(StaticMatrix<typename TestFixture::Kernel, 0, 2>(0.));
		problem[0].b(StaticVector<typename TestFixture::Kernel, 0>(0.));

		// TODO: can we make lbx() etc. accept initializer lists?
		problem[0].lbx(typename TestFixture::Vector {-10000., -inf<Real>(), -10000.});
		problem[0].ubx(typename TestFixture::Vector {10000., inf<Real>(), 10000.});
		problem[0].lbu(typename TestFixture::Vector {-inf<Real>(), -10000.});
		problem[0].ubu(typename TestFixture::Vector {inf<Real>(), 10000.});

		ws.solve();
		auto solution = ws.solution();

		using Vector = typename TestFixture::Vector;

		EXPECT_PRED2(ApproxEqual(1e-6), solution[0].x(), (Vector {
			146.566682434017,
			-427.218558989821,
			-345.347969700289
		}));

		EXPECT_PRED2(ApproxEqual(1e-6), solution[0].u(), (Vector {
			769.663140469139,
			-191.122524763114
		}));

		EXPECT_PRED2(ApproxEqual(1e-6), solution[1].x(), (Vector {}));
		EXPECT_PRED2(ApproxEqual(1e-6), solution[1].u(), (Vector {}));
	}

	TYPED_TEST_P(TreeQpWorkspaceSolveTest, testSolve2)
	{
		using Real = typename TestFixture::Real;
		using Kernel = typename TestFixture::Kernel;
		
		typename TestFixture::Workspace workspace {OcpSize {3, 0, 0}, OcpSize {0, 0, 0}};
		
		auto& stage0 = workspace.problem()[0];
		stage0.gaussNewtonCostApproximation(
			Vector {1., 2., 42.},
			IdentityMatrix<Kernel> {3u},
			Matrix {3u, 0u}
		);
		stage0.stateBounds(-inf<Real>(), inf<Real>());
		stage0.inputBounds(-inf<Real>(), inf<Real>());

		workspace.solve();
	
		EXPECT_EQ(forcePrint(workspace.solution()[0].x()), forcePrint(Vector {-1., -2., -42.}));
	}
#endif

	REGISTER_TYPED_TEST_SUITE_P(TreeQpWorkspaceSolveTest,
		//testMoveConstructor, 
		testSolve0
		, testSolve0Unconstrained
		, testSolveDimitrisRandom
		//,testSolve1
		//DISABLED_testSolve1stage1d, 
		//testSolve2stage1d, 
		//testSolve2stage5d, 
		//testInfiniteBoundsAll, 
		//testInfiniteBoundsSome,
		//testSolve2
	);
}