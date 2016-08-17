#include <integrator/rk4.hpp>
#include <casadi_interface/GeneratedFunction.hpp>

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

class PendulumODE
{
public:
	static unsigned const NX = 2;
	static unsigned const NU = 1;

	typedef Eigen::Matrix<double, NX, 1> StateVector;
	typedef Eigen::Matrix<double, NU, 1> InputVector;
	typedef Eigen::Matrix<double, NX, NX, Eigen::ColMajor> StateStateMatrix;
	typedef Eigen::Matrix<double, NX, NU, Eigen::ColMajor> StateInputMatrix;

	void operator()(double t, StateVector const& x0, InputVector const& u0,	StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B) const
	{
		static casadi_interface::GeneratedFunction const _ode(CASADI_GENERATED_FUNCTION_INTERFACE(pendulum_ode));
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), nullptr, A.data(), B.data(), nullptr, nullptr});
	}

	void operator()(double t, StateVector const& x0, InputVector const& u0,
			StateVector const& x0_seed, InputVector const& u_seed, StateVector& xdot, StateVector& xdot_sens) const
	{
		static casadi_interface::GeneratedFunction const _ode(CASADI_GENERATED_FUNCTION_INTERFACE(pendulum_ode_sens));
		_ode({&t, x0.data(), u0.data(), x0_seed.data(), u_seed.data()}, {xdot.data(), xdot_sens.data()});
	}

	StateVector operator()(double t, StateVector const& x0, InputVector const& u0) const
	{
		static casadi_interface::GeneratedFunction const _ode(CASADI_GENERATED_FUNCTION_INTERFACE(pendulum_ode));

		StateVector xdot;
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), nullptr, nullptr, nullptr, nullptr, nullptr});

		return xdot;
	}
};

class rk4_test : public ::testing::Test
{
protected:
	typedef PendulumODE ODE;
	typedef tmpc::RK4 Integrator;

	ODE ode_;
	Integrator integrator_ {0.01};

	std::ifstream test_data_ {"test/data/rk4/pendulum.txt"};

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
		ODE::StateVector x0;
		ODE::InputVector u;
		ODE::StateVector xplus;
		double qf;
		ODE::StateStateMatrix A;
		ODE::StateInputMatrix B;
		ODE::StateVector qA;
		ODE::InputVector qB;

		friend std::istream& operator>>(std::istream& is, TestPoint& p)
		{
			return is >> p.t >> p.x0 >> p.u >> p.xdot >> p.Aode >> p.Bode >> p.xplus
						>> p.qf >> p.A >> p.B >> p.qA >> p.qB;
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

		EXPECT_THAT(as_container(xplus), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.xplus)));
		EXPECT_THAT(as_container(A), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.A)));
		EXPECT_THAT(as_container(B), testing::Pointwise(FloatNearPointwise(1e-5), as_container(p.B)));

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
