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
	typedef Eigen::Matrix<double, NX + NU, 1> StateInputVector;
	typedef Eigen::Matrix<double, NX, NX + NU, Eigen::ColMajor> ODEJacobianMatrix;

	PendulumODE() : _ode(CASADI_GENERATED_FUNCTION_INTERFACE(pendulum_ode)) {}

	void ODE(double t, StateInputVector const& z0, StateVector& xdot, ODEJacobianMatrix& J) const
	{
		_ode({&t, z0.data()}, {xdot.data(), J.data()});
	}

private:
	casadi_interface::GeneratedFunction const _ode;
};

TEST(rk4_test, integrate)
{
	typedef PendulumODE ODE;
	typedef camels::RK4<PendulumODE> Integrator;

	ODE ode;
	Integrator integrator(ode, 0.01);

	std::ifstream test_data("data/rk4/pendulum.txt");

	double t;
	ODE::StateVector xdot_expected;
	ODE::ODEJacobianMatrix Jode_expected;
	Integrator::StateInputVector z0;
	Integrator::StateVector xplus_expected;
	Integrator::ODEJacobianMatrix J_expected;

	unsigned count = 0;
	while (test_data >> t >> z0 >> xdot_expected >> Jode_expected >> xplus_expected >> J_expected)
	{
		{
			ODE::StateVector xdot;
			ODE::ODEJacobianMatrix Jode;
			ode.ODE(t, z0, xdot, Jode);

			EXPECT_TRUE(xdot.isApprox(xdot_expected));
			EXPECT_TRUE(Jode.isApprox(Jode_expected));
		}

		{
			Integrator::StateVector xplus;
			Integrator::ODEJacobianMatrix J;
			integrator.Integrate(t, z0, xplus, J);

			EXPECT_TRUE(xplus.isApprox(xplus_expected));
			EXPECT_TRUE(J    .isApprox(J_expected    ));
		}

		++count;
	}

	EXPECT_EQ(count, 600);
}
