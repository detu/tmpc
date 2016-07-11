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
	typedef Eigen::Matrix<double, NX + NU, 1> StateInputVector;
	typedef Eigen::Matrix<double, NX, NX + NU, Eigen::ColMajor> ODEJacobianMatrix;
	typedef Eigen::Matrix<double, NX, NX, Eigen::ColMajor> StateStateMatrix;
	typedef Eigen::Matrix<double, NX, NU, Eigen::ColMajor> StateInputMatrix;

	PendulumODE() {}

	void ODE(double t, StateInputVector const& z0, StateVector& xdot, ODEJacobianMatrix& J) const
	{
		static casadi_interface::GeneratedFunction const _ode(CASADI_GENERATED_FUNCTION_INTERFACE(pendulum_ode_jac));
		_ode({&t, z0.data()}, {xdot.data(), J.data()});
	}

	void ODE(double t, StateVector const& x0, InputVector const& u0, StateVector& xdot, StateStateMatrix& A, StateInputMatrix& B) const
	{
		static casadi_interface::GeneratedFunction const _ode(CASADI_GENERATED_FUNCTION_INTERFACE(pendulum_ode_AB));
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data()});
	}
};

class rk4_test : public ::testing::Test
{
protected:
	typedef PendulumODE ODE;
	typedef tmpc::RK4<PendulumODE> Integrator;

	ODE ode_;
	Integrator integrator_ {ode_, 0.01};

	std::ifstream test_data_ {"data/rk4/pendulum.txt"};
};

TEST_F(rk4_test, integrate_works)
{
	double t;
	ODE::StateVector xdot_expected;
	ODE::ODEJacobianMatrix Jode_expected;
	Integrator::StateInputVector z0;
	Integrator::StateVector xplus_expected;
	Integrator::ODEJacobianMatrix J_expected;

	unsigned count = 0;
	while (test_data_ >> t >> z0 >> xdot_expected >> Jode_expected >> xplus_expected >> J_expected)
	{
		{
			ODE::StateVector xdot;
			ODE::ODEJacobianMatrix Jode;
			ode_.ODE(t, z0, xdot, Jode);

			EXPECT_TRUE(xdot.isApprox(xdot_expected));
			EXPECT_TRUE(Jode.isApprox(Jode_expected));
		}

		{
			Integrator::StateVector xplus;
			Integrator::ODEJacobianMatrix J;
			integrator_.Integrate(t, z0, xplus, J);

			EXPECT_TRUE(xplus.isApprox(xplus_expected));
			EXPECT_TRUE(J    .isApprox(J_expected    ));
		}

		++count;
	}

	EXPECT_EQ(count, 600);
}

TEST_F(rk4_test, integrate_new_interface_works)
{
	double t;
	ODE::StateVector xdot_expected;
	ODE::ODEJacobianMatrix Jode_expected;
	Integrator::StateInputVector z0;
	Integrator::StateVector xplus_expected;
	Integrator::ODEJacobianMatrix J_expected;

	unsigned count = 0;
	while (test_data_ >> t >> z0 >> xdot_expected >> Jode_expected >> xplus_expected >> J_expected)
	{
		{
			ODE::StateVector xdot;
			ODE::StateStateMatrix A;
			ODE::StateInputMatrix B;
			ode_.ODE(t, tmpc::top_rows<ODE::NX>(z0), tmpc::bottom_rows<ODE::NU>(z0), xdot, A, B);

			EXPECT_TRUE(xdot.isApprox(xdot_expected));
			EXPECT_TRUE(A.isApprox(tmpc::left_cols <ODE::NX>(Jode_expected)));
			EXPECT_TRUE(B.isApprox(tmpc::right_cols<ODE::NU>(Jode_expected)));
		}

		{
			Integrator::StateVector xplus;
			ODE::StateStateMatrix A;
			ODE::StateInputMatrix B;
			integrator_.Integrate(t, tmpc::top_rows<ODE::NX>(z0), tmpc::bottom_rows<ODE::NU>(z0), xplus, A, B);

			EXPECT_TRUE(xplus.isApprox(xplus_expected));
			EXPECT_TRUE(A.isApprox(tmpc::left_cols <ODE::NX>(J_expected)));
			EXPECT_TRUE(B.isApprox(tmpc::right_cols<ODE::NU>(J_expected)));
		}

		++count;
	}

	EXPECT_EQ(count, 600);
}
