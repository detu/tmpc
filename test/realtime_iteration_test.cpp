#include "../include/core/RealtimeIteration.hpp"
#include "../include/qp/CondensingSolver.hpp"
#include "../include/integrator/RK4.hpp"
#include "../include/core/OptimalControlProblem.hpp"

#include <Eigen/Dense>

#include <gtest/gtest.h>

struct Dimensions
{
	static unsigned const NX = 2;
	static unsigned const NU = 1;
	static unsigned const ND = 0;
	static unsigned const NDT = 0;
};

class SampleOCP : public tmpc::OptimalControlProblem<SampleOCP, Dimensions::NX, Dimensions::NU>
{
public:
	static unsigned const NX  = Dimensions::NX;
	static unsigned const NU  = Dimensions::NU;
	//static unsigned const NZ
	static unsigned const NC  = 0;
	static unsigned const NCT = 0;

	typedef double Scalar;
	typedef Eigen::Matrix<Scalar, NC, 1> ConstraintVector;
	typedef Eigen::Matrix<Scalar, NC, NW> ConstraintJacobianMatrix;
	typedef Eigen::Matrix<Scalar, NCT, 1> TerminalConstraintVector;
	typedef Eigen::Matrix<Scalar, NCT, NX> TerminalConstraintJacobianMatrix;

	SampleOCP(std::size_t nt)
	:	tmpc::OptimalControlProblem<SampleOCP, Dimensions::NX, Dimensions::NU>(nt)
	{
		_x_min << -1., -1.;
		_x_max <<  1.,  1.;
		_u_min << -1.;
		_u_max <<  1.;
		_x_terminal_min << -1., -1.;
		_x_terminal_max <<  1.,  1.;
	}

	void LagrangeTerm(unsigned i, StateInputVector const& z, StateInputVector& g, LagrangeHessianMatrix& H) const
	{
		H << 66,  78,  90,
			 78,  93, 108,
			 90, 108, 126;

		g << 0, 0, 0;
	}

	StateVector const& getStateMin() const { return _x_min; }
	StateVector const& getStateMax() const { return _x_max; }
	StateVector const& getTerminalStateMin() const { return _x_terminal_min; }
	StateVector const& getTerminalStateMax() const { return _x_terminal_max; }

	InputVector const& getInputMin() const { return _u_min;	}
	InputVector const& getInputMax() const { return _u_max;	}

	void MayerTerm(const StateVector& x, StateVector& g, MayerHessianMatrix& H) const
	{
		H << 10, 14,
		     14, 20;

		g << 0, 0;
	}

	void PathConstraints(unsigned i, const StateInputVector& z,
		ConstraintJacobianMatrix& D, ConstraintVector& d_min, ConstraintVector& d_max) const
	{
	}

	void TerminalConstraints(const StateVector& x, TerminalConstraintJacobianMatrix& D,
		TerminalConstraintVector& d_min, TerminalConstraintVector& d_max) const
	{
	}

private:
	StateVector _x_min;
	StateVector _x_max;
	StateVector _x_terminal_min;
	StateVector _x_terminal_max;
	InputVector _u_min;
	InputVector _u_max;
};

typedef SampleOCP OCP;

class DiscreteTimeModel
{
public:
	typedef OCP::StateVector StateVector;
	typedef OCP::StateInputVector StateInputVector;
	typedef OCP::ODEJacobianMatrix ODEJacobianMatrix;

	static unsigned const NX = OCP::NX;
	static unsigned const NU = OCP::NU;

	DiscreteTimeModel() {}

	void Integrate(double t0, StateInputVector const& z0, StateVector& x_next, ODEJacobianMatrix& J) const
	{
		J << 1.,  1., 0.5,
		     0.,  1., 1. ;

		x_next = J * z0;
	}

	double timeStep() const { return 1.; }
};

class Controller
{
	typedef DiscreteTimeModel Integrator;
	typedef tmpc::CondensingSolver<OCP::NX, OCP::NU, OCP::NC, OCP::NCT> QPSolver;
	typedef tmpc::RealtimeIteration<OCP, Integrator, QPSolver> RealtimeIteration;
	typedef RealtimeIteration::Trajectory Trajectory;

	OCP _ocp;
	Integrator _integrator;
	QPSolver _qpSolver;
	RealtimeIteration _rti;

public:
	Controller(unsigned Nt = 2)
	:	_ocp(Nt)
	,	_qpSolver(Nt, camels::qpOASESOptions::MPC().setPrintLevel(qpOASES::PL_LOW))
	,	_rti(_ocp, _integrator, _qpSolver,
			tmpc::ConstantTrajectory<Trajectory>(Nt, OCP::StateVector::Zero(), OCP::InputVector::Zero()))
	{
	}

	void Preparation()
	{
		_rti.Preparation();
	}

	OCP::InputVector Feedback(OCP::StateVector const& x0)
	{
		return _rti.Feedback(x0);
	}
};

TEST(realtime_iteration, sample_ocp)
{
	Controller controller;

	OCP::StateVector x0;
	x0 << 1, 0;

	controller.Preparation();
	auto const u0 = controller.Feedback(x0);

	OCP::InputVector u0_expected;
	u0_expected << -0.690877362606266;
	EXPECT_TRUE(u0.isApprox(u0_expected, 1e-6));

	std::cout << u0;
}
