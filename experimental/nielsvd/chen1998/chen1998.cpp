#include "core/RealtimeIteration.hpp"
#include "core/gauss_newton.hpp"
#include "qp/CondensingSolver.hpp"
#include "qp/Condensing.hpp"
#include "qp/qpOASESProgram.hpp"
#include "integrator/RK4.hpp"
#include "casadi_interface/GeneratedFunction.hpp"

#include <Eigen/Dense>

#include "chen1998_ode_generated.h"

#include <iostream>

using namespace tmpc;

class ODE
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
		static casadi_interface::GeneratedFunction const _ode(CASADI_GENERATED_FUNCTION_INTERFACE(chen1998_ode_AB));
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), A.data(), B.data()});
	}

	StateVector operator()(double t, StateVector const& x0, InputVector const& u0) const
	{
		static casadi_interface::GeneratedFunction const _ode(CASADI_GENERATED_FUNCTION_INTERFACE(chen1998_ode_AB));

		StateVector xdot;
		_ode({&t, x0.data(), u0.data()}, {xdot.data(), nullptr, nullptr});

		return xdot;
	}
};

class SampleOCP
{
public:
	static unsigned const NX  = 2;
	static unsigned const NU  = 1;
	static unsigned const NC  = 0;
	static unsigned const NCT = 0;

	typedef Eigen::Matrix<double, NX, 1> StateVector;
	typedef Eigen::Matrix<double, NU, 1> InputVector;

	SampleOCP(std::size_t nt)
	:	nt_(nt)
	{
		_x_min << -100., -100.;
		_x_max <<  100.,  100.;
		_u_min << -2.;
		_u_max <<  2.;
		_x_terminal_min << -100., -100.;
		_x_terminal_max <<  100.,  100.;
	}

	std::size_t getNumberOfIntervals() const
	{
		return nt_;
	}

	ODE const& getODE() const
	{
		return _ode;
	}

	template <typename StateVector, typename InputVector, typename QMatrix, typename RMatrix, typename SMatrix,
		typename StateGradientVector, typename InputGradientVector>
	void LagrangeTerm(unsigned i, StateVector const& x, InputVector const& u, QMatrix& Q, RMatrix& R, SMatrix& S,
			StateGradientVector& q, InputGradientVector& r) const
	{
    Eigen::Matrix<double, NX+NU, 1> res;
    Eigen::Matrix<double, NX+NU, NX+NU> W;
    Eigen::Matrix<double, NX+NU, NX> C;
    Eigen::Matrix<double, NX+NU, NU> D;
    res << x,
           u;
    W << 0.5,  0, 0,
			   0,  0.5, 0,
         0,    0, 1;
    C << 1,    0,
         0,    1,
         0,    0;
    D << 0,
         0,
         1;
    Gauss_Newton_approximation(res, W, C, D, Q, R, S, q, r);
	}

	StateVector const& getStateMin() const { return _x_min; }
	StateVector const& getStateMax() const { return _x_max; }
	StateVector const& getTerminalStateMin() const { return _x_terminal_min; }
	StateVector const& getTerminalStateMax() const { return _x_terminal_max; }

	InputVector const& getInputMin() const { return _u_min;	}
	InputVector const& getInputMax() const { return _u_max;	}

	template <typename StateVector, typename GradientVector, typename HessianMatrix>
	void MayerTerm(const StateVector& x, GradientVector& g, HessianMatrix& H) const
	{
		H << 16.5926, 11.5926,
		     11.5926, 16.5926;

		g = transpose(x)*H;
	}

	template <typename StateInputVector, typename ConstraintJacobianMatrix, typename ConstraintVector>
	void PathConstraints(unsigned i, const StateInputVector& z,
		ConstraintJacobianMatrix& D, ConstraintVector& d_min, ConstraintVector& d_max) const
	{
	}

	template <typename StateVector, typename TerminalConstraintJacobianMatrix, typename TerminalConstraintVector>
	void TerminalConstraints(const StateVector& x, TerminalConstraintJacobianMatrix& D,
		TerminalConstraintVector& d_min, TerminalConstraintVector& d_max) const
	{
	}

private:
	std::size_t nt_;
	StateVector _x_min;
	StateVector _x_max;
	StateVector _x_terminal_min;
	StateVector _x_terminal_max;
	InputVector _u_min;
	InputVector _u_max;
  ODE _ode;
};

typedef SampleOCP OCP;

int main() {
  // Consensing QP solver
  typedef CondensingSolver<OCP::NX, OCP::NU, OCP::NC, OCP::NCT> QPSolver;
  // Integrator
  typedef RK4 Integrator;
  // Real Time Iteration
  typedef RealtimeIteration<OCP, Integrator, QPSolver> RealtimeIteration;
  typedef typename RealtimeIteration::WorkingPoint WorkingPoint;

  const int N = 13;
  const double dt = 0.1;

  QPSolver qp(N);
  OCP ocp(N);
  Integrator integrator(dt);
  RealtimeIteration rti(ocp,integrator,qp,WorkingPoint(N, OCP::StateVector::Zero(), OCP::InputVector::Zero()));

  OCP::StateVector x;
  OCP::InputVector u;

  // x << 0.292, 0.228;
  x << -0.683, -0.864;
  for(int i=0; i<10; ++i) {
    u = rti.Feedback(x);
    rti.Preparation();
  }

  std::cout << rti.getWorkingPoint() << std::endl;

  std::cout << "Einde.\n";
}
