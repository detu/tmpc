#include <integrator/RK4.hpp>
#include <casadi_interface/GeneratedFunction.hpp>

#include "pendulum_ode_generated.h"

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include <fstream>

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
	typedef Eigen::Matrix<double, 0 , 1> ParamVector;
	typedef Eigen::Matrix<double, NX, NX, Eigen::ColMajor> StateStateMatrix;
	typedef Eigen::Matrix<double, NX, NU, Eigen::ColMajor> StateInputMatrix;

	void operator()(double t, StateVector const& x0, InputVector const& u0,
			ParamVector const&, StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B) const
	{
		static casadi_interface::GeneratedFunction const _ode(CASADI_GENERATED_FUNCTION_INTERFACE(pendulum_ode_AB));
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data()});
	}
};

class rk4_test : public ::testing::Test
{
protected:
	typedef PendulumODE ODE;
	typedef tmpc::RK4 Integrator;

	ODE ode_;
	Integrator integrator_ {0.01};

	std::ifstream test_data_ {"data/rk4/pendulum.txt"};
};

TEST_F(rk4_test, integrate_new_interface_works)
{
	double t;
	ODE::StateVector xdot_expected;
	ODE::StateStateMatrix Aode_expected;
	ODE::StateInputMatrix Bode_expected;
	ODE::StateVector x0;
	ODE::InputVector u;
	ODE::StateVector xplus_expected;
	ODE::StateStateMatrix A_expected;
	ODE::StateInputMatrix B_expected;
	ODE::ParamVector p;

	unsigned count = 0;
	while (test_data_ >> t >> x0 >> u >> xdot_expected >> Aode_expected >> Bode_expected >> xplus_expected >> A_expected >> B_expected)
	{
		{
			ODE::StateVector xdot;
			ODE::StateStateMatrix A;
			ODE::StateInputMatrix B;
			ode_(t, x0, u, p, xdot, A, B);

			EXPECT_TRUE(xdot.isApprox(xdot_expected));
			EXPECT_TRUE(A.isApprox(Aode_expected));
			EXPECT_TRUE(B.isApprox(Bode_expected));
		}

		{
			ODE::StateVector xplus;
			ODE::StateStateMatrix A;
			ODE::StateInputMatrix B;
			integrator_.Integrate(ode_, t, x0, u, p, xplus, A, B);

			EXPECT_TRUE(xplus.isApprox(xplus_expected));
			EXPECT_TRUE(A.isApprox(A_expected));
			EXPECT_TRUE(B.isApprox(B_expected));
		}

		++count;
	}

	EXPECT_EQ(count, 600);
}
