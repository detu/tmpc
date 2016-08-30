#include "../include/core/RealtimeIteration.hpp"
#include "../include/qp/CondensingSolver.hpp"
#include "../include/qp/HPMPCSolver.hpp"
#include "../include/integrator/rk4.hpp"

#include <Eigen/Dense>

#include <gtest/gtest.h>

#include <type_traits>

struct Dimensions
{
	static unsigned const NX = 2;
	static unsigned const NU = 1;
	static unsigned const ND = 0;
	static unsigned const NDT = 0;
};

// An empty class for ODE, since we provide the discrete-time dynamics ourselves.
struct ODE
{
	ODE() {}
};

class SampleOCP
{
public:
	static unsigned const NX  = Dimensions::NX;
	static unsigned const NU  = Dimensions::NU;
	static unsigned const NC  = 0;
	static unsigned const NCT = 0;

	typedef Eigen::Matrix<double, NX, 1> StateVector;
	typedef Eigen::Matrix<double, NU, 1> InputVector;

	SampleOCP(std::size_t nt)
	:	nt_(nt)
	{
		_x_min << -1., -1.;
		_x_max <<  1.,  1.;
		_u_min << -1.;
		_u_max <<  1.;
		_x_terminal_min << -1., -1.;
		_x_terminal_max <<  1.,  1.;
	}

	std::size_t getNumberOfIntervals() const
	{
		return nt_;
	}

	ODE const& getODE() const
	{
		static ODE const ode;
		return ode;
	}

	template <typename StateVector, typename InputVector, typename QMatrix, typename RMatrix, typename SMatrix,
		typename StateGradientVector, typename InputGradientVector>
	void LagrangeTerm(unsigned i, StateVector const& x, InputVector const& u, QMatrix& Q, RMatrix& R, SMatrix& S,
			StateGradientVector& q, InputGradientVector& r) const
	{
		Q << 66,  78,
			 78,  93;
		R << 126;
		S << 90,
			108;
		q << 0., 0.;
		r << 0.;
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
		H << 10, 14,
		     14, 20;

		g << 0, 0;
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
};

typedef SampleOCP OCP;

class DiscreteTimeModel
{
public:
	double timeStep() const { return 1.; }
};

template <typename StateVector, typename InputVector, typename NextStateVector, typename AMatrix, typename BMatrix>
void integrate(DiscreteTimeModel const& integrator, ODE const&, double t0, StateVector const& x0, InputVector const& u, NextStateVector& x_next,
		Eigen::MatrixBase<AMatrix>& A, Eigen::MatrixBase<BMatrix>& B)
{
	A << 1.,  1.,
		 0.,  1.;

	B << 0.5,
		 1. ;

	x_next = A * x0 + B * u;
}

template <typename RealtimeIteration>
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
	typedef typename RealtimeIteration::Integrator Integrator;
	//typedef tmpc::RealtimeIteration<OCP, Integrator, QPSolver> RealtimeIteration;
	typedef typename RealtimeIteration::WorkingPoint WorkingPoint;
	typedef typename RealtimeIteration::QPSolver QPSolver;

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
		tmpc::RealtimeIteration<OCP, DiscreteTimeModel, tmpc::CondensingSolver>
,		tmpc::RealtimeIteration<OCP, DiscreteTimeModel, tmpc::HPMPCSolver     >
	> RTITypes;

TYPED_TEST_CASE(RealtimeIterationTest, RTITypes);

TYPED_TEST(RealtimeIterationTest, GivesCorrectU0)
{
	OCP::StateVector x;
	OCP::InputVector u;

	// Step 0
	{
		x << 1, 0;
		u = this->Feedback(x);
		this->Preparation();

		OCP::InputVector u_expected;
		u_expected << -0.690877362606266;

		EXPECT_TRUE(u.isApprox(u_expected, 1e-6));
	}

	// Step 1
	{
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
	}
}
