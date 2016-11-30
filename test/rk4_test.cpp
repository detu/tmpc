#include <tmpc/integrator/rk4.hpp>
#include <tmpc/casadi_interface/GeneratedFunction.hpp>

#include "pendulum_ode_generated.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Eigen/Dense>

#include <fstream>
#include "gtest_tools_eigen.hpp"

template<typename Matrix>
std::istream& operator>>(std::istream& is, Eigen::MatrixBase<Matrix>& m)
{
	for (typename Matrix::Index i = 0; i < m.rows(); ++i)
		for (typename Matrix::Index j = 0; j < m.cols(); ++j)
			is >> m(i, j);

	return is;
}

class PendulumODEBase
{
public:
	static unsigned const NX = 2;
	static unsigned const NU = 1;
	static unsigned const NQ = 2;
	static unsigned const NR = 2;

	typedef Eigen::Matrix<double, NX, 1> StateVector;
	typedef Eigen::Matrix<double, NU, 1> InputVector;
	typedef Eigen::Matrix<double, NQ, 1> QuadVector;
	typedef Eigen::Matrix<double, NR, 1> ResVector;
	typedef Eigen::Matrix<double, NX, NX, Eigen::ColMajor> StateStateMatrix;
	typedef Eigen::Matrix<double, NX, NU, Eigen::ColMajor> StateInputMatrix;
	typedef Eigen::Matrix<double, NU, NU, Eigen::ColMajor> InputInputMatrix;
	typedef Eigen::Matrix<double, NQ, NX, Eigen::ColMajor> QuadStateMatrix;
	typedef Eigen::Matrix<double, NQ, NU, Eigen::ColMajor> QuadInputMatrix;
	typedef Eigen::Matrix<double, NR, NX, Eigen::ColMajor> ResStateMatrix;
	typedef Eigen::Matrix<double, NR, NU, Eigen::ColMajor> ResInputMatrix;

protected:
	casadi_interface::GeneratedFunction<CASADI_GENERATED_FUNCTION_INTERFACE(pendulum_ode), 3, 9> const _ode;
};

class PendulumODE : public PendulumODEBase
{
public:
	/**
	 * \brief Evaluates ODE.
	 */
	void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B) const
	{
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data(), nullptr, nullptr, nullptr, nullptr, nullptr, nullptr});
	}

	/**
	 * \brief Evaluates ODE and quadrature.
	 */
	void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
		QuadVector& q, QuadStateMatrix& qA, QuadInputMatrix& qB) const
	{
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data(), q.data(), qA.data(), qB.data(), nullptr, nullptr, nullptr});
	}

	/**
	 * \brief Evaluates ODE, quadrature and residuals.
	 */
	void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
		QuadVector& q, QuadStateMatrix& qA, QuadInputMatrix& qB, ResVector& r, ResStateMatrix& rA, ResInputMatrix& rB) const
	{
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data(), q.data(), qA.data(), qB.data(), r.data(), rA.data(), rB.data()});
	}

	void operator()(double t, StateVector const& x0, InputVector const& u0,
			StateVector const& x0_seed, InputVector const& u_seed, StateVector& xdot, StateVector& xdot_sens) const
	{
		static casadi_interface::GeneratedFunction<CASADI_GENERATED_FUNCTION_INTERFACE(pendulum_ode_sens), 5, 2> const _ode;
		_ode({&t, x0.data(), u0.data(), x0_seed.data(), u_seed.data()}, {xdot.data(), xdot_sens.data()});
	}

	/**
	 * \brief Evaluates ODE without sensitivities.
	 */
	StateVector operator()(double t, StateVector const& x0, InputVector const& u0) const
	{
		StateVector xdot;
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr});

		return xdot;
	}
};

class PendulumODE_r : public PendulumODEBase
{
public:
	/**
	 * \brief Evaluates ODE and residuals.
	 */
	void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B,
		ResVector& r, ResStateMatrix& rA, ResInputMatrix& rB) const
	{
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data(), nullptr, nullptr, nullptr, r.data(), rA.data(), rB.data()});
	}
};

class rk4_test : public ::testing::Test
{
protected:
	typedef PendulumODEBase ODE;
	typedef tmpc::RK4 Integrator;

	PendulumODE ode_;
	PendulumODE_r ode_r_;
	Integrator integrator_ {0.01};

	std::ifstream test_data_ {std::string(TEST_DATA_PATH) + "/rk4/pendulum.txt"};

	void SetUp() override
	{
		ASSERT_TRUE(test_data_);
	}

	struct TestPoint
	{
		double t;
		ODE::StateVector xdot;
		ODE::StateStateMatrix Aode;
		ODE::StateInputMatrix Bode;
		ODE::QuadVector q;
		ODE::QuadStateMatrix qA_ode;
		ODE::QuadInputMatrix qB_ode;
		ODE::StateVector x0;
		ODE::InputVector u;
		ODE::StateVector xplus;
		ODE::QuadVector qf;
		ODE::StateStateMatrix A;
		ODE::StateInputMatrix B;
		ODE::QuadStateMatrix qA;
		ODE::QuadInputMatrix qB;
		ODE::ResVector r;
		ODE::ResStateMatrix rA_ode;
		ODE::ResInputMatrix rB_ode;
		double cf;
		ODE::StateVector cA;
		ODE::InputVector cB;
		ODE::StateStateMatrix cQ;
		ODE::InputInputMatrix cR;
		ODE::StateInputMatrix cS;

		friend std::istream& operator>>(std::istream& is, TestPoint& p)
		{
			// keys = ['t', 'x0', 'u', 'xdot', 'A_ode', 'B_ode', 'q', 'qA_ode', 'qB_ode', 'x_plus', 'A', 'B', 'qf', 'qA', 'qB']
			return is >> p.t >> p.x0 >> p.u >> p.xdot  >> p.Aode   >> p.Bode
					                        >> p.q     >> p.qA_ode >> p.qB_ode
					                        >> p.xplus >> p.A      >> p.B
						                    >> p.qf    >> p.qA     >> p.qB
											>> p.r	   >> p.rA_ode >> p.rB_ode
											>> p.cf    >> p.cA     >> p.cB
											>> p.cQ    >> p.cR     >> p.cS;
		};
	};
};

TEST_F(rk4_test, ode_correct)
{
	TestPoint p;

	unsigned count = 0;
	while (test_data_ >> p)
	{
		ODE::StateVector xdot;
		ODE::StateStateMatrix A;
		ODE::StateInputMatrix B;
		ode_(p.t, p.x0, p.u, xdot, A, B);

		EXPECT_TRUE(xdot.isApprox(p.xdot));
		EXPECT_TRUE(A.isApprox(p.Aode));
		EXPECT_TRUE(B.isApprox(p.Bode));

		++count;
	}

	EXPECT_EQ(count, 600);
}

TEST_F(rk4_test, ode_q_correct)
{
	TestPoint p;

	unsigned count = 0;
	while (test_data_ >> p)
	{
		ODE::StateVector xdot;
		ODE::StateStateMatrix A;
		ODE::StateInputMatrix B;
		ODE::QuadVector q;
		ODE::QuadStateMatrix qA;
		ODE::QuadInputMatrix qB;
		ode_(p.t, p.x0, p.u, xdot, A, B, q, qA, qB);

		EXPECT_TRUE(xdot.isApprox(p.xdot));
		EXPECT_TRUE(A.isApprox(p.Aode));
		EXPECT_TRUE(B.isApprox(p.Bode));
		EXPECT_TRUE(q.isApprox(p.q));
		EXPECT_TRUE(qA.isApprox(p.qA_ode));
		EXPECT_TRUE(qB.isApprox(p.qB_ode));

		++count;
	}

	EXPECT_EQ(count, 600);
}

TEST_F(rk4_test, ode_qr_correct)
{
	TestPoint p;

	unsigned count = 0;
	while (test_data_ >> p)
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

		EXPECT_TRUE(xdot.isApprox(p.xdot));
		EXPECT_TRUE(A.isApprox(p.Aode));
		EXPECT_TRUE(B.isApprox(p.Bode));
		EXPECT_TRUE(q.isApprox(p.q));
		EXPECT_TRUE(qA.isApprox(p.qA_ode));
		EXPECT_TRUE(qB.isApprox(p.qB_ode));
		EXPECT_TRUE(r.isApprox(p.r));
		EXPECT_TRUE(rA.isApprox(p.rA_ode));
		EXPECT_TRUE(rB.isApprox(p.rB_ode));

		++count;
	}

	EXPECT_EQ(count, 600);
}

TEST_F(rk4_test, integrate_correct)
{
	TestPoint p;

	unsigned count = 0;
	while (test_data_ >> p)
	{
		ODE::StateVector xplus;
		ODE::StateStateMatrix A;
		ODE::StateInputMatrix B;
		integrate(integrator_, ode_, p.t, p.x0, p.u, xplus, A, B);

		/*
		EXPECT_EQ(print_wrap(xplus), print_wrap(p.xplus));
		EXPECT_EQ(print_wrap(A), print_wrap(p.A));
		EXPECT_EQ(print_wrap(B), print_wrap(p.B));
		*/
		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.xplus)));
		EXPECT_THAT(as_container(A), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.A)));
		EXPECT_THAT(as_container(B), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.B)));

		++count;
	}

	EXPECT_EQ(count, 600);
}

TEST_F(rk4_test, integrate_q_correct)
{
	TestPoint p;

	unsigned count = 0;
	while (test_data_ >> p)
	{
		ODE::StateVector xplus;
		ODE::StateStateMatrix A;
		ODE::StateInputMatrix B;
		ODE::QuadVector qf;
		ODE::QuadStateMatrix qA;
		ODE::QuadInputMatrix qB;
		integrate(integrator_, ode_, p.t, p.x0, p.u, xplus, A, B, qf, qA, qB);

		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));
		EXPECT_THAT(as_container(A    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.A    )));
		EXPECT_THAT(as_container(B    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.B    )));
		EXPECT_THAT(as_container(qf   ), testing::Pointwise(FloatNearPointwise(1e-4), as_container(p.qf   )));
		EXPECT_THAT(as_container(qA   ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.qA   )));
		EXPECT_THAT(as_container(qB   ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.qB   )));

		++count;
	}

	EXPECT_EQ(count, 600);
}

TEST_F(rk4_test, integrate_r_correct)
{
	TestPoint p;

	unsigned count = 0;
	while (test_data_ >> p)
	{
		ODE::StateVector xplus;
		ODE::StateStateMatrix A;
		ODE::StateInputMatrix B;
		double cf;
		ODE::StateVector cA;
		ODE::InputVector cB;
		ODE::StateStateMatrix cQ;
		ODE::InputInputMatrix cR;
		ODE::StateInputMatrix cS;
		integrate(integrator_, ode_r_, p.t, p.x0, p.u, xplus, A, B, cf, cA, cB, cQ, cR, cS);

		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));
		EXPECT_THAT(as_container(A    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.A    )));
		EXPECT_THAT(as_container(B    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.B    )));
		EXPECT_NEAR(cf, p.cf, 1e-10);
		EXPECT_THAT(as_container(cA   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cA   )));
		EXPECT_THAT(as_container(cB   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cB   )));
		EXPECT_THAT(as_container(cQ   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cQ   )));
		EXPECT_THAT(as_container(cR   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cR   )));
		EXPECT_THAT(as_container(cS   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cS   )));
		++count;
	}

	EXPECT_EQ(count, 600);
}

TEST_F(rk4_test, integrate_qr_correct)
{
	TestPoint p;

	unsigned count = 0;
	while (test_data_ >> p)
	{
		ODE::StateVector xplus;
		ODE::StateStateMatrix A;
		ODE::StateInputMatrix B;
		ODE::QuadVector qf;
		ODE::QuadStateMatrix qA;
		ODE::QuadInputMatrix qB;
		double cf;
		ODE::StateVector cA;
		ODE::InputVector cB;
		ODE::StateStateMatrix cQ;
		ODE::InputInputMatrix cR;
		ODE::StateInputMatrix cS;
		integrate(integrator_, ode_, p.t, p.x0, p.u, xplus, A, B, qf, qA, qB, cf, cA, cB, cQ, cR, cS);

		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));
		EXPECT_THAT(as_container(A    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.A    )));
		EXPECT_THAT(as_container(B    ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.B    )));
		EXPECT_THAT(as_container(qf   ), testing::Pointwise(FloatNearPointwise(1e-4), as_container(p.qf   )));
		EXPECT_THAT(as_container(qA   ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.qA   )));
		EXPECT_THAT(as_container(qB   ), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.qB   )));
		EXPECT_NEAR(cf, p.cf, 1e-10);
		EXPECT_THAT(as_container(cA   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cA   )));
		EXPECT_THAT(as_container(cB   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cB   )));
		EXPECT_THAT(as_container(cQ   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cQ   )));
		EXPECT_THAT(as_container(cR   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cR   )));
		EXPECT_THAT(as_container(cS   ), testing::Pointwise(FloatNearPointwise(1e-10), as_container(p.cS   )));
		++count;
	}

	EXPECT_EQ(count, 600);
}

TEST_F(rk4_test, integrate_no_sens_correct)
{
	TestPoint p;

	unsigned count = 0;
	while (test_data_ >> p)
	{
		auto const xplus = integrate(integrator_, ode_, p.t, p.x0, p.u);
		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));

		++count;
	}

	EXPECT_EQ(count, 600);
}

TEST_F(rk4_test, integrate_fd_correct)
{
	TestPoint p;

	unsigned count = 0;
	while (test_data_ >> p)
	{
		ODE::StateVector xplus;
		ODE::StateStateMatrix A;
		ODE::StateInputMatrix B;
		integrate(integrator_, ode_, p.t, p.x0, p.u, xplus, A, B);

		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));
		EXPECT_THAT(as_container(A), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.A)));
		EXPECT_THAT(as_container(B), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.B)));

		++count;
	}

	EXPECT_EQ(count, 600);
}
