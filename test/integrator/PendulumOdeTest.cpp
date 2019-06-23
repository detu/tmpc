#include "PendulumOde.hpp"
#include "PendulumData.hpp"

#include <tmpc/Testing.hpp>


namespace tmpc :: testing
{
	class PendulumOdeTest 
	: 	public Test
	{
	protected:
		using Real = double;
		typedef PendulumOdeBase ODE;

		PendulumOde ode_;

		std::vector<TestPoint> test_data_;

		void SetUp() override
		{
			test_data_ = loadPendulumData();
			ASSERT_EQ(test_data_.size(), 600);
		}
	};


	TEST_F(PendulumOdeTest, testOdeCorrect)
	{
		for (auto p : test_data_)
		{
			ODE::StateVector xdot;
			ODE::StateStateMatrix A;
			ODE::StateInputMatrix B;
			ode_(p.t, p.x0, p.u, xdot, A, B);

			double const tol = 1e-6;
			EXPECT_TRUE(approxEqual(xdot, p.xdot, tol));
			EXPECT_TRUE(approxEqual(A, p.Aode, tol));
			EXPECT_TRUE(approxEqual(B, p.Bode, tol));
		}
	}
	

	TEST_F(PendulumOdeTest, testOdeQCorrect)
	{
		for (auto p : test_data_)
		{
			ODE::StateVector xdot;
			ODE::StateStateMatrix A;
			ODE::StateInputMatrix B;
			ODE::QuadVector q;
			ODE::QuadStateMatrix qA;
			ODE::QuadInputMatrix qB;
			ode_(p.t, p.x0, p.u, xdot, A, B, q, qA, qB);

			double const tol = 1e-6;
			EXPECT_TRUE(approxEqual(xdot, p.xdot, tol));
			EXPECT_TRUE(approxEqual(A, p.Aode, tol));
			EXPECT_TRUE(approxEqual(B, p.Bode, tol));
			EXPECT_TRUE(approxEqual(q, p.q, tol));
			EXPECT_TRUE(approxEqual(qA, p.qA_ode, tol));
			EXPECT_TRUE(approxEqual(qB, p.qB_ode, tol));
		}
	}


	TEST_F(PendulumOdeTest, testOdeQrCorrect)
	{
		for (auto p : test_data_)
		{
			ODE::StateVector xdot;
			ODE::StateStateMatrix A;
			ODE::StateInputMatrix B;
			ODE::QuadVector q;
			ODE::QuadStateMatrix qA;
			ODE::QuadInputMatrix qB;
			ODE::ResVector r;
			ODE::ResStateMatrix rA;
			ODE::ResInputMatrix rB;
			ode_(p.t, p.x0, p.u, xdot, A, B, q, qA, qB, r, rA, rB);

			double const tol = 1e-6;
			EXPECT_TRUE(approxEqual(xdot, p.xdot, tol));
			EXPECT_TRUE(approxEqual(A, p.Aode, tol));
			EXPECT_TRUE(approxEqual(B, p.Bode, tol));
			EXPECT_TRUE(approxEqual(q, p.q, tol));
			EXPECT_TRUE(approxEqual(qA, p.qA_ode, tol));
			EXPECT_TRUE(approxEqual(qB, p.qB_ode, tol));
			EXPECT_TRUE(approxEqual(r, p.r, tol));
			EXPECT_TRUE(approxEqual(rA, p.rA_ode, tol));
			EXPECT_TRUE(approxEqual(rB, p.rB_ode, tol));
		}
	}
}