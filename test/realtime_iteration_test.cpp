#include "../include/core/RealtimeIteration.hpp"
#include "../include/qp/CondensingSolver.hpp"
#include "../include/qp/HPMPCSolver.hpp"
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

template <typename QPSolver>
class RealtimeIterationTest : public ::testing::Test
{
public:
	RealtimeIterationTest(unsigned Nt = 2)
	:	_ocp(Nt)
	,	_qpSolver(Nt)
	,	_rti(_ocp, _integrator, _qpSolver, WorkingPoint(Nt, OCP::StateVector::Zero(), OCP::InputVector::Zero()))
	{
	}

protected:
	typedef DiscreteTimeModel Integrator;
	//typedef tmpc::CondensingSolver<OCP::NX, OCP::NU, OCP::NC, OCP::NCT> QPSolver;
	typedef tmpc::RealtimeIteration<OCP, Integrator, QPSolver> RealtimeIteration;
	typedef typename RealtimeIteration::WorkingPoint WorkingPoint;

	OCP _ocp;
	Integrator _integrator;
	QPSolver _qpSolver;
	RealtimeIteration _rti;

	void Preparation()
	{
		_rti.Preparation();
	}

	OCP::InputVector Feedback(OCP::StateVector const& x0)
	{
		return _rti.Feedback(x0);
	}
};

typedef ::testing::Types<
		tmpc::CondensingSolver<OCP::NX, OCP::NU, OCP::NC, OCP::NCT>
,		tmpc::HPMPCSolver     <OCP::NX, OCP::NU, OCP::NC, OCP::NCT>
	> QPSolvers;

TYPED_TEST_CASE(RealtimeIterationTest, QPSolvers);

TYPED_TEST(RealtimeIterationTest, GivesCorrectU0)
{
	OCP::StateVector x;
	OCP::InputVector u;

	// Step 0
	{
		this->Preparation();

		x << 1, 0;
		u = this->Feedback(x);

		OCP::InputVector u_expected;
		u_expected << -0.690877362606266;
		EXPECT_TRUE(u.isApprox(u_expected, 1e-6));

		std::cout << "u[0] = " << u.transpose() << std::endl;
	}

	// Step 1
	{
		this->Preparation();

		std::cout << "Working point before step 1: " << std::endl;
		std::cout << this->_rti.getWorkingPoint() << std::endl;

		{
			OCP::StateVector x;
			OCP::InputVector u;

			x << 0.654561318696867,	 -0.690877362606266;	u << 0.215679569867116;
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_x(0).isApprox(x, 1e-6));
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_u(0).isApprox(u, 1e-6));

			x << 0.0715237410241597, -0.475197792739149;	u << 0.215679569867116;
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_x(1).isApprox(x, 1e-6));
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_u(1).isApprox(u, 1e-6));

			x << 0.0715237410241597, -0.475197792739149;
			EXPECT_TRUE(this->_rti.getWorkingPoint().get_x(2).isApprox(x, 1e-6));
		}

		x << 0.654561318696867,	-0.690877362606266;
		u = this->Feedback(x);

		OCP::InputVector u_expected;
		u_expected << 0.218183;
		EXPECT_TRUE(u.isApprox(u_expected, 1e-5));

		std::cout << "u[1] = " << u.transpose() << std::endl;
	}
}
