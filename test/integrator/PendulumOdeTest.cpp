#include "PendulumOde.hpp"
#include "PendulumData.hpp"

#include <tmpc/test_tools.hpp>


namespace tmpc :: testing
{
	class PendulumOdeTest 
	: 	public ::testing::Test
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

			MatrixApproxEquality const is_approx(1e-6);
			EXPECT_PRED2(is_approx, xdot, p.xdot);
			EXPECT_PRED2(is_approx, A, p.Aode);
			EXPECT_PRED2(is_approx, B, p.Bode);
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

			MatrixApproxEquality const is_approx(1e-6);
			EXPECT_PRED2(is_approx, xdot, p.xdot);
			EXPECT_PRED2(is_approx, A, p.Aode);
			EXPECT_PRED2(is_approx, B, p.Bode);
			EXPECT_PRED2(is_approx, q, p.q);
			EXPECT_PRED2(is_approx, qA, p.qA_ode);
			EXPECT_PRED2(is_approx, qB, p.qB_ode);
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

			MatrixApproxEquality const is_approx(1e-6);
			EXPECT_PRED2(is_approx, xdot, p.xdot);
			EXPECT_PRED2(is_approx, A, p.Aode);
			EXPECT_PRED2(is_approx, B, p.Bode);
			EXPECT_PRED2(is_approx, q, p.q);
			EXPECT_PRED2(is_approx, qA, p.qA_ode);
			EXPECT_PRED2(is_approx, qB, p.qB_ode);
			EXPECT_PRED2(is_approx, r, p.r);
			EXPECT_PRED2(is_approx, rA, p.rA_ode);
			EXPECT_PRED2(is_approx, rB, p.rB_ode);
		}
	}
}